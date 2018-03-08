function uploadFunction() {
    var fileInput = $('#files');
    if (!window.FileReader) {
        alert('Your browser is not supported')
    }
    var input = fileInput.get(0);

    // Create a reader object
    var reader = new FileReader();
    if (input.files.length) {
        var textFile = input.files[0];
        reader.readAsText(textFile);
        $(reader).on('load', processFile);
    } else {
        alert('Please upload a file before continuing')
    }
}

function processFile(e) {
    var file = e.target.result,
        results;
    if (file && file.length) {
        $('#lab').text(file);
    }
}

function myFunction() {
    alert("I am an alert box!");
}
