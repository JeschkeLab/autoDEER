from PyQt6.QtWidgets import QApplication
from autodeer import autoDEERUI

app = QApplication([])
window = autoDEERUI()
window.show()
app.exec()