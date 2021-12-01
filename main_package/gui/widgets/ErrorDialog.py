from PyQt5.QtWidgets import QMessageBox

# Have to call ErrorDialog.show() to show it.
class ErrorDialog(QMessageBox):
    def __init__(self):
        super().__init__(icon=QMessageBox.Warning)

    def show_with_message(self, message: str) -> None:
        bulletpoint: str = b'\xe2\x80\xa2'.decode()
        rightarrow: str = b'\xe2\x9e\xa4'.decode()
        delta: str = b'\xce\x94'.decode()
        lam: str = b'\xce\xbb'.decode()

        message = f'Error:\n{message}\n\n' \
                  f'Possible causes:\n\n' \
                  f'{bulletpoint} Integration step size {delta}t is too large:\n' \
                  f'{rightarrow} Choose a smaller {delta}t and restart the simulation.\n\n' \
                  f'{bulletpoint} One of the {lam}\'s is too large:\n' \
                  f'{rightarrow} Choose a smaller {lam} and restart the simulation.\n\n' \
                  f'{bulletpoint} One of the mSq\'s is negative:\n' \
                  f'{rightarrow} Choose a positive mSq and restart the simulation.\n'

        self.setText(message)
        self.show()