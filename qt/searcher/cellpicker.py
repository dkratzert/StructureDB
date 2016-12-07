'''
Created on 09.02.2015

 ----------------------------------------------------------------------------
* "THE BEER-WARE LICENSE" (Revision 42):
* <daniel.kratzert@uni-freiburg.de> wrote this file. As long as you retain this 
* notice you can do whatever you want with this stuff. If we meet some day, and 
* you think this stuff is worth it, you can buy me a beer in return.
* ----------------------------------------------------------------------------

@author: daniel
'''


def get_cif_cell(filename):
    '''
    Get the cell of a cif file
    _cell_length_a                    29.227(4)
    _cell_length_b                    6.6568(8)
    _cell_length_c                    11.8204(14)
    _cell_angle_alpha                 90
    _cell_angle_beta                  107.055(5)
    _cell_angle_gamma                 90
    ...
    _shelx_hkl_file
    ;
    hkl data
    ;
    '''
    cells = []
    name, a, b, c, alpha, beta, gamma = None, None, None, None, None, None, None
    with open(filename) as f:
        cell = [None, None, None, None, None, None, None]
        for line in f:
            if line.startswith('data_'):
                name = line.split('_')[1].strip('\n')
                cell[0] = name
            if line.startswith('_cell_length_a'):
                a = line.split()[1].split('(')[0]
                cell[1] = a
            if line.startswith('_cell_length_b'):
                b = line.split()[1].split('(')[0]
                cell[2] = b
            if line.startswith('_cell_length_c'):
                c = line.split()[1].split('(')[0]
                cell[3] = c
            if line.startswith('_cell_angle_alpha'):
                alpha = line.split()[1].split('(')[0]
                cell[4] = alpha
            if line.startswith('_cell_angle_beta'):
                beta = line.split()[1].split('(')[0]
                cell[5] = beta
            if line.startswith('_cell_angle_gamma'):
                gamma = line.split()[1].split('(')[0]
                cell[6] = gamma
            if all(cell[n] for n, i in enumerate(cell)):
            #if cell[0] and cell[1] and cell[2] and cell[3] and cell[4] and cell[5] and cell[6]:
                cells.append(cell[:])
                cell = [None, None, None, None, None, None, None]               
    return cells

    

def get_res_cell(filename):
    '''
    Returns the unit cell parameters from the list file as list:
    ['a', 'b', 'c', 'alpha', 'beta', 'gamma']
    '''
    from misc import open_file_read
    file_list = open_file_read(filename)
    cell = False
    for line in file_list:
        if line.startswith('CELL'):
            cell = line.split()[2:8]
            try:
                cell = [float(i) for i in cell]
                if not len(cell) == 6:
                    raise ValueError
            except(ValueError):
                print('Bad cell parameters in {0}.'.format(filename))
                return False
            #break
    if not cell:
        #print('Unable to find unit cell parameters in the file.')
        return False
    return cell


if __name__ == '__main__':
    cells = get_cif_cell('c:/temp/c2c_final.cif')
    cellres = get_res_cell('c:/temp/c2c.res')
    for i in cells:
        print(i)
    print(cellres)