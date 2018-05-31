
scale_ratio = 1.2

class Position():
    def __init__(self, x, y):
        self.x = x
        self.y = y

class Rectangle():
    def __init__(self, x, y, dx, dy, tooltip_url, css_class='rectangle'):
        self.pos = Position(x * 1000, (y) * scale_ratio)
        self.width = max(1, dx * 1000 - 2)
        self.height = dy - 1
        self.tooltip_url = tooltip_url
        self.css_class = css_class

class Arrow():
    def __init__(self, left_x, y, right_direction=True, css_class='arrow'):
        self.left_x = left_x * 1000
        self.y = y * scale_ratio
        self.right_direction = right_direction
        self.css_class = css_class

class HorizontalLine():
    def __init__(self, left_x, right_x, y, height=2, css_style='horizontal_line'):
        self.css_class = css_style
        self.left_x = left_x * 1000
        self.width = (right_x - left_x) * 1000
        self.y = (y) * scale_ratio
        self.height = height

class VerticalLine():
    def __init__(self, x, top_y, bottom_y, width=2, css_style='vertical_line'):
        self.css_class = css_style
        self.x = x * 1000
        self.top_y = (top_y) * scale_ratio
        self.height = (bottom_y - top_y) * scale_ratio
        self.width = width
        if css_style == 'vertical_dotted_line' or css_style == 'vertical_dashed_line':
            self.width = 0

class Text():
    def __init__(self, text, left_x, top_y, css_class='text', width=80, height=80):
        self.css_class = css_class
        self.text = text
        self.left_x = left_x * 1000
        self.top_y = top_y * scale_ratio
        self.width = 20 + (len(str(text)) - 1) * 7.5
        self.height = 25 + (len(str(text)) - 1) * 7.5

class CycleSection():
    def __init__(self, top_y, bottom_y, cycle_name, rectangles, segments_labels, background_opacity, css_class=''):
        self.css_class = css_class
        self.rectangles = rectangles
        self.segments_labels = segments_labels
        top_y *= scale_ratio
        bottom_y *= scale_ratio
        self.top_y = top_y + 5
        self.name = cycle_name
        self.height = bottom_y - top_y + 5
        self.center_y = (bottom_y + top_y) / 2 - 2
        self.background_opacity = background_opacity
