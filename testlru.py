from pylru import lrudecorator


@lrudecorator(4)
def moo(label, val):
    print "label"
    return val*val

print moo("10",10)
print moo("10",10)
print moo("10",10)
print moo("10",10)
print moo("10",10)
print moo("11",11)
print moo("12",12)
print moo("13",13)
print moo("10",10)
print moo("14",14)
print moo("15",15)
print moo("10",10)
print moo("11",11)
