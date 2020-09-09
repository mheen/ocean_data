from datetime import datetime
import traceback

class Type:
    def __init__(self, prefix):
        self.prefix = prefix

ERROR = Type("ERROR")
INFO = Type("")

def info(log_file, message):
    _write(log_file, INFO, message)

def error(log_file, message, exception=None):
    _write(log_file, ERROR, message)

def _write(log_file, type, message, ex=None):
    timestamp = datetime.now().strftime('%d-%m-%Y %H:%M:%S')
    log_msg = f'{timestamp} - {type.prefix}: {message}'
    if ex:
        log_msg += '\n' + str(traceback.format_exception(etype=type(ex), value=ex, tb=ex.__traceback__))
    print(log_msg)
    if log_file is not None:
        with open(log_file,'a') as f:
            f.write(log_msg + '\n')
        