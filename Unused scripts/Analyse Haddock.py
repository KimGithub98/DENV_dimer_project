import sys
print(sys.argv)

pdbs = sys.argv[1:]
Para_struct = '~/software/haddock2.4-2021-01/run_' + pdb + '/structures/it1/water/analysis/'
True_struct = '~/software/haddock2.4-2021-01/run_' + pdb + '_T/structures/it1/water/analysis/'
haddock_matrix = [["pdb", "Evdw", "Eelec", "Edesolv", "H_score", "hbond"]]

if len(sys.argv) > 2:
    for pdb in pdbs:
        for struct in [Para_struct, True_struct]:
            energies = open(struct + 'energies.disp')
            edesolv = open(struct + 'edesolv.disp')
            hbond = open(struct + 'ana_hbonds.lis')
            #pdb = open(struct + 'haddock2.4-2021-01_ava.pdb')
            Evdw = float(energies.readlines()[175][8:17])
            Eelec = float(energies.readlines()[176][8:17])
            Edesolv = float(energies.readlines()[55][12:21])
            H_score = 1.0*Evdw + 0.2*Eelec + 1.0 *Edesolv
            hbond_list = []
            for x in hbond.read().split("\n"):
                hbond_i = x[48:50] + x[53:56] + "-" + x[74:76] + x[79:82]
                hbond_list.append(hbond_i)
            if struct == Para_struct:
                haddock_matrix.append([pdb + "_P", Evdw, Eelec, Edesolv, H_score, hbond_list])
            elif struct == True_struct:
                haddock_matrix.append([pdb + "_T", Evdw, Eelec, Edesolv, H_score, hbond_list])
print(haddock_matrix)
