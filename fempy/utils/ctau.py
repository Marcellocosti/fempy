c = 299792458 # m/s#
hbarc = 197.3269804 # MeV·fm

def get_ctau_time(time):
    return time*c

def get_ctau_width(width):
    return hbarc / width

mean_lives = { # s
    411: 1033e-15,
    421: 410.3e-15,
}

widths = { # MeV
    413: 0.0834,
    3224: 36.2
}

print(get_ctau_time(mean_lives[411]))
print(get_ctau_time(mean_lives[421]))
print(get_ctau_width(widths[413]))
