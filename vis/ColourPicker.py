from PySide2.QtWidgets import *
from PySide2.QtGui import *
from PySide2.QtCore import *

class ColourPicker():
    kelly_colour_vals = (0x222222, 0xf3c300, 0x875692, 0xf38400, 0xa1caf1, 0xbe0032, 0xc2b280, 0x848482, 0x008856, 0xe68fac, 0x0067a5, 0xf99379, 0x604e97, 0xf6a600, 0xb3446c, 0xdcd300, 0x882d17, 0x8db600, 0x654522, 0xe25822, 0x2b3d26) #note: white has been removed

    colours = [QColor.fromRgb(c) for c in kelly_colour_vals]

    def __init__(self):
        self.next_index = 0
        self.cache = {} #key -> colour_index

    def get(self, key):
        if key in self.cache:
            index = self.cache[key]
            colour = ColourPicker.colours[index]
        else:
            colour = ColourPicker.colours[self.next_index]
            self.cache[key] = self.next_index
            self.next_index = (self.next_index + 1) % len(ColourPicker.colours)

        return colour

    @Slot()
    def reset(self):
        self.cache.clear()
        self.next_index = 0
            
            
