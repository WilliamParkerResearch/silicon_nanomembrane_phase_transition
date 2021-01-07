atomic_calculations = {
    1: {"pseudopotential": "Si.pz-n-rrkjus_psl.0.1.UPF",        "total energy": -10.98894304},
    2: {"pseudopotential": "Si.pbe-nl-kjpaw_psl.1.0.0.UPF",     "total energy": -46.39167521},
    3: {"pseudopotential": "Si.pbesol-nl-kjpaw_psl.1.0.0.UPF",  "total energy": -45.27537292},
    4: {"pseudopotential": "Si.SCAN.UPF",                       "total energy":  -7.53623525}
}

for id, info in atomic_calculations.items():
    print("\nCalculation ID:", id)
    for key in info:
        if key == "total energy":
            print(key + ':', info[key], 'Ry')
        else:
            print(key + ':', info[key])
