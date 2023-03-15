# Interface to MOSEK for AmpliconArchitect for Python3
#
# Supports all versions of MOSEK >= 8
#
import logging
import sys

import mosek

# Check Mosek version
mosek_ver = mosek.Env.getversion()
logging.info("Mosek version is {}".format('.'.join([str(x) for x in mosek_ver])))
mosek_major = mosek_ver[0]

if sys.version_info < (3, 0) and mosek_major >= 10:
    logging.warning("Mosek version is " + '.'.join([str(x) for x in mosek_ver]) +
                    " which requires python3. Exiting.\n")
    sys.exit(1)


# MOSEK logging
mosek_logger = logging.getLogger('MOSEK')


def moseklogfunc(msg):
    mosek_logger.debug(msg.rstrip())


class fusionlogger:
    def write(self, msg):
        moseklogfunc(msg)

    def flush(self):
        pass


# Calls MOSEK to solve one instance of the problem
def call_mosek(n, m, asub, aval, coeff_c, coeff_f, coeff_g, const_h):
    mosek_logger.info("Beginning MOSEK call")

    ## Enable this line to ALWAYS save all Mosek inputs
    #save_mosek_input(n, m, asub, aval, coeff_c, coeff_f, coeff_g, const_h)

    try:
        # Determine which MOSEK routing to call
        if mosek_major == 8:
            return call_mosek_scopt(n, m, asub, aval, coeff_c, coeff_f, coeff_g, const_h)
        elif mosek_major == 9:
            return call_mosek_fusion(n, m, asub, aval, coeff_c, coeff_f)
        elif mosek_major >= 10:
            return call_mosek_acc(n, m, asub, aval, coeff_c, coeff_f)
        else:
            raise Exception("Unsupported MOSEK version {}".format(mosek_major))
    except Exception as e:
        # If an error occurred in the MOSEK call then save
        # all input data to a JSON file so they can be loaded
        # to recreate the MOSEK problem in a stand-alone way.
        mosek_logger.error("Error when using MOSEK: {}".format(e))
        print("Error when using MOSEK: {}".format(e))
        filename = save_mosek_input(n, m, asub, aval, coeff_c, coeff_f, coeff_g, const_h)        
        mosek_logger.info("Saved MOSEK inputs to {}. Submit that file to support to reproduce the issue.".format(filename))
        raise e


'''
This method works with MOSEK == 8.

Solves the problem

minimize     c^T * x - sum_i(f_i * log(g_i * x_i + h_i))
subject to   A * x == 0
             x >= 0
'''             
def call_mosek_scopt(n, m, asub, aval, coeff_c, coeff_f, coeff_g, const_h):
    with mosek.Env() as env:
        with env.Task() as task:
            task.set_Stream(mosek.streamtype.log, moseklogfunc)

            numvar = n + m
            numcon = 2 * n 
            
            task.appendcons(numcon)
            task.appendvars(numvar)

            task.putvarboundslice(0, numvar, [mosek.boundkey.lo]* numvar, [0.0]*numvar, [0.0]*numvar)
            task.putconboundslice(0, numcon, [mosek.boundkey.fx]* numcon, [0.0]*numcon, [0.0]*numcon)

            for i in range(numcon):
                task.putarow(i, asub[i], aval[i])

            task.putclist(range(numvar), coeff_c)

            task.putobjsense(mosek.objsense.minimize)

            task.putSCeval([mosek.scopr.log] * (n + m), range(n + m), coeff_f, coeff_g, const_h)

            task.optimize()
            task.solutionsummary(mosek.streamtype.log)

            if task.getsolsta(mosek.soltype.itr) != mosek.solsta.optimal:
                raise Exception("Failed to solve to optimality. Solution status {}".format(task.getsolsta(mosek.soltype.itr)))

            res = [ 0.0 ] * numvar
            task.getsolutionslice(mosek.soltype.itr, mosek.solitem.xx, 0, numvar, res)

            return res


'''
This method works with MOSEK >= 10.

Solves the problem

minimize     c^T * x - sum_i(f_i * log(x_i))
subject to   A * x == 0

Comments, compared to MOSEK 8 model:

We ignore the normalizing coefficient h_i from log(g_i * x_i + h_i) and consider only log(g_i * x_i). 
Subject to that change we can also skip g_i since it only changes the constant term in the objective.
The condition x>=0 is implicit by x appearing in the logarithm.
'''             
def call_mosek_acc(n, m, asub, aval, coeff_c, coeff_f):
    with mosek.Task() as task:
        task.set_Stream(mosek.streamtype.log, moseklogfunc)

        task.appendvars(2*(n + m))
        task.appendcons(2*n)
        task.putvarboundsliceconst(0, 2*(n + m), mosek.boundkey.fr, 0, 0)

        for i in range(2*n):
            task.putarow(i, asub[i], aval[i])

        task.putconboundsliceconst(0, 2*n, mosek.boundkey.fx, 0, 0)

        task.appendafes(2*(n + m) + 1)
        task.putafefentrylist(range(0, 2*(n + m)), range(0, 2*(n + m)), [1.0]*(2*(n + m)))
        task.putafeg(2*(n + m), 1.0)

        expdom = task.appendprimalexpconedomain()
        task.appendaccs([expdom]*(n + m), sum([[i, 2*(n + m), i + n + m] for i in range(n + m)], []), None)

        task.putclist(range(0, n + m), coeff_c)
        task.putclist(range(n + m, 2 * (n + m)), coeff_f)
        
        task.putobjsense(mosek.objsense.minimize)

        task.optimize()
        task.solutionsummary(mosek.streamtype.log)
        
        if task.getsolsta(mosek.soltype.itr) != mosek.solsta.optimal:
            raise Exception("Failed to solve to optimality. Solution status {}".format(task.getsolsta(mosek.soltype.itr)))

        return task.getxxslice(mosek.soltype.itr, 0, n + m)


'''
This method works with MOSEK >= 9.

Solves the problem

minimize     c^T * x - sum_i(f_i * log(x_i))
subject to   A * x == 0

Comments, compared to MOSEK 10 model:

A simple model in the higher level MOSEK Fusion. Anyhow, we do not expect MOSEK 9 users, really. 
Either stay with MOSEK 8 or otherwise there is no reason not to upgrade all the way to MOSEK 10. 

This model can be used in MOSEK >= 9, but it invokes the additional Fusion modeling layer,
which the model from call_mosek_acc skips. If it behaves well though, we could make
it the default. It should be fast enough, and is more readable.
'''             
def call_mosek_fusion(n, m, asub, aval, coeff_c, coeff_f):
    from mosek.fusion import Model, Domain, Expr, Matrix, ObjectiveSense, SolutionStatus
    with Model() as M:
        M.setLogHandler(fusionlogger())

        x = M.variable(n + m)
        t = M.variable(n + m)

        for i in range(2*n):
            M.constraint(Expr.dot(aval[i], x.pick(asub[i])), Domain.equalsTo(0))

        M.constraint(Expr.hstack(x, Expr.constTerm(n + m, 1.0), t), Domain.inPExpCone())

        M.objective(ObjectiveSense.Minimize, Expr.add(Expr.dot(coeff_c, x), Expr.dot(coeff_f, t)))

        M.solve()

        if M.getPrimalSolutionStatus() != SolutionStatus.Optimal:
            raise Exception("Failed to solve to optimality. Solution status {}".format(M.getPrimalSolutionStatus()))

        return x.level()


# Debug functions. Dumping input data.
mosek_save_num = 1
def save_mosek_input(n, m, asub, aval, coeff_c, coeff_f, coeff_g, const_h):
    import json
    global mosek_save_num
    filename = "mosekinput-{}.json".format(mosek_save_num)
    data = { "n": n, "m": m, 
             "asub": asub, "aval": aval, 
             "coeff_c": coeff_c, 
             "coeff_f": coeff_f, "coeff_g": coeff_g, "const_h": const_h}
    
    with open(filename, "w") as f:
        json.dump(data, f)
    
    mosek_save_num += 1
    return filename


# Debug functions. Loading input data.
def load_mosek_input(filename):
    import json
    with open(filename, "r") as f:
        data = json.load(f)
    return data["n"], data["m"], data["asub"], data["aval"], data["coeff_c"], data["coeff_f"], data["coeff_g"], data["const_h"]
