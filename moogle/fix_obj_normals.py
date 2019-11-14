
base = '/home/rnakanishi/git/bahamut-lib/results/particles/3d/enright/'
folders = [
    # 'pls128_1',
    'semiL128_2', 'weno128_2', 'pls128_2'
]

for folder in folders:
    for f in range(180, 181):
        fpath = base + folder + '/mesh/{:04d}.obj'.format(f)
        correct = open(base + folder + '/mesh/n{:04d}.obj'.format(f), 'w')
        with open(fpath) as fp:
            line = fp.readline()
            while line:
                if(line[0] == 'f'):
                    line = line[2:-2]
                    # values = [int(val) for val in line.split(' ')]
                    newline = "f "
                    newline += " ".join(["%s//%s" % (value, value)
                                         for value in line.split(' ')])
                    line = newline + '\n'
                correct.write(line)
                line = fp.readline()
        print("Read file " + folder + " %04d" % f)
        correct.close()
