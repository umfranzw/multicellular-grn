from PySide2.QtWidgets import QApplication
from MainWindow import MainWindow

def main():
    app = QApplication([])
    win = MainWindow()
    win.show()
    app.exec_()

main()
