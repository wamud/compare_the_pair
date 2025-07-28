PLUS = lambda x, y: x + y
MINUS = lambda x, y: x - y


class Err(float):
    def __init__(self, err): 
        self.val = val
        self.err = err
        self.sign = PLUS 

    def __add__(self, other):
        return self.sign() super().__add__(self.val, other)



def max_err(fn, *errs): 
    def _wrap(*args):
        args = list(args)
        for i, err in enumerate(err): 
            if err is not None:
                arg[i] = err(args[i]) 
        return fn(*args) 


