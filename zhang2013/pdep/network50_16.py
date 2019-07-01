species(
    label = 'C[CH]CCCC(49)',
    structure = SMILES('C[CH]CCCC'),
    E0 = (2.58175,'kJ/mol'),
    modes = [
        HarmonicOscillator(frequencies=([2750,2770,2790,2810,2830,2850,1350,1400,1450,1500,700,800,1000,1100,1350,1400,900,1100,2750,2770,2790,2810,2830,2850,1425,1437.5,1450,1225,1250,1275,1270,1305,1340,700,750,800,300,350,400,3025,407.5,1350,352.5,200,800,1600],'cm^-1')),
        HinderedRotor(inertia=(0.156089,'amu*angstrom^2'), symmetry=1, barrier=(3.5888,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.156089,'amu*angstrom^2'), symmetry=1, barrier=(3.5888,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.156089,'amu*angstrom^2'), symmetry=1, barrier=(3.5888,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.156089,'amu*angstrom^2'), symmetry=1, barrier=(3.5888,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.156089,'amu*angstrom^2'), symmetry=1, barrier=(3.5888,'kJ/mol'), semiclassical=False),
    ],
    spinMultiplicity = 2,
    opticalIsomers = 1,
    molecularWeight = (85.1674,'amu'),
    collisionModel = TransportData(shapeIndex=2, epsilon=(3202.63,'J/mol'), sigma=(6.0508,'angstroms'), dipoleMoment=(0,'C*m'), polarizability=(0,'angstroms^3'), rotrelaxcollnum=0, comment="""Epsilon & sigma estimated with Tc=500.24 K, Pc=32.8 bar (from Joback method)"""),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[1.10039,0.0591158,-2.83098e-05,5.96613e-09,-4.29309e-13,418.934,27.2711], Tmin=(100,'K'), Tmax=(1919.43,'K')), NASAPolynomial(coeffs=[15.9939,0.0337697,-1.29499e-05,2.17603e-09,-1.36862e-13,-6346.85,-57.0233], Tmin=(1919.43,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(2.58175,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(436.51,'J/(mol*K)'), comment="""Thermo group additivity estimation: group(Cs-CsCsHH) + group(Cs-CsCsHH) + group(Cs-CsCsHH) + group(Cs-CsCsHH) + group(Cs-CsHHH) + group(Cs-CsHHH) + radical(RCCJC)"""),
)

species(
    label = '[CH2]CCCCC(10)',
    structure = SMILES('[CH2]CCCCC'),
    E0 = (13.3818,'kJ/mol'),
    modes = [
        HarmonicOscillator(frequencies=([2750,2800,2850,1350,1500,750,1050,1375,1000,2750,2764.29,2778.57,2792.86,2807.14,2821.43,2835.71,2850,1425,1433.33,1441.67,1450,1225,1241.67,1258.33,1275,1270,1293.33,1316.67,1340,700,733.333,766.667,800,300,333.333,366.667,400,3000,3100,440,815,1455,1000,200,800,1600],'cm^-1')),
        HinderedRotor(inertia=(0.156089,'amu*angstrom^2'), symmetry=1, barrier=(3.5888,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.156089,'amu*angstrom^2'), symmetry=1, barrier=(3.5888,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.156089,'amu*angstrom^2'), symmetry=1, barrier=(3.5888,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.156089,'amu*angstrom^2'), symmetry=1, barrier=(3.5888,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.156089,'amu*angstrom^2'), symmetry=1, barrier=(3.5888,'kJ/mol'), semiclassical=False),
    ],
    spinMultiplicity = 2,
    opticalIsomers = 1,
    molecularWeight = (85.1674,'amu'),
    collisionModel = TransportData(shapeIndex=2, epsilon=(3202.63,'J/mol'), sigma=(6.0508,'angstroms'), dipoleMoment=(0,'C*m'), polarizability=(0,'angstroms^3'), rotrelaxcollnum=0, comment="""Epsilon & sigma estimated with Tc=500.24 K, Pc=32.8 bar (from Joback method)"""),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[0.795713,0.061583,-2.70843e-05,7.8915e-10,1.81285e-12,1732.05,27.3958], Tmin=(100,'K'), Tmax=(1257.3,'K')), NASAPolynomial(coeffs=[11.7154,0.0395171,-1.58794e-05,2.86532e-09,-1.94139e-13,-2015.57,-31.7693], Tmin=(1257.3,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(13.3818,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(436.51,'J/(mol*K)'), comment="""Thermo group additivity estimation: group(Cs-CsCsHH) + group(Cs-CsCsHH) + group(Cs-CsCsHH) + group(Cs-CsCsHH) + group(Cs-CsHHH) + group(Cs-CsHHH) + radical(RCCJ)"""),
)

species(
    label = 'CC[CH]CCC(50)',
    structure = SMILES('CC[CH]CCC'),
    E0 = (2.59378,'kJ/mol'),
    modes = [
        HarmonicOscillator(frequencies=([2750,2770,2790,2810,2830,2850,1350,1400,1450,1500,700,800,1000,1100,1350,1400,900,1100,2750,2770,2790,2810,2830,2850,1425,1437.5,1450,1225,1250,1275,1270,1305,1340,700,750,800,300,350,400,3025,407.5,1350,352.5,200.549,803.991,1611.32],'cm^-1')),
        HinderedRotor(inertia=(0.153001,'amu*angstrom^2'), symmetry=1, barrier=(3.55299,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.153001,'amu*angstrom^2'), symmetry=1, barrier=(3.55299,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.153001,'amu*angstrom^2'), symmetry=1, barrier=(3.55299,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.153001,'amu*angstrom^2'), symmetry=1, barrier=(3.55299,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.153001,'amu*angstrom^2'), symmetry=1, barrier=(3.55299,'kJ/mol'), semiclassical=False),
    ],
    spinMultiplicity = 2,
    opticalIsomers = 1,
    molecularWeight = (85.1674,'amu'),
    collisionModel = TransportData(shapeIndex=2, epsilon=(3202.63,'J/mol'), sigma=(6.0508,'angstroms'), dipoleMoment=(0,'C*m'), polarizability=(0,'angstroms^3'), rotrelaxcollnum=0, comment="""Epsilon & sigma estimated with Tc=500.24 K, Pc=32.8 bar (from Joback method)"""),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[1.33686,0.0573244,-2.58738e-05,4.74284e-09,-2.3902e-13,408.839,26.2442], Tmin=(100,'K'), Tmax=(2003.81,'K')), NASAPolynomial(coeffs=[17.0745,0.0329272,-1.28644e-05,2.16251e-09,-1.35163e-13,-7307.24,-64.1353], Tmin=(2003.81,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(2.59378,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(436.51,'J/(mol*K)'), comment="""Thermo group additivity estimation: group(Cs-CsCsHH) + group(Cs-CsCsHH) + group(Cs-CsCsHH) + group(Cs-CsCsHH) + group(Cs-CsHHH) + group(Cs-CsHHH) + radical(RCCJCC)"""),
)

species(
    label = '[CH2]C(C)CCC(68)',
    structure = SMILES('[CH2]C(C)CCC'),
    E0 = (7.5235,'kJ/mol'),
    modes = [
        HarmonicOscillator(frequencies=([1380,1390,370,380,2900,435,2750,2783.33,2816.67,2850,1425,1450,1225,1275,1270,1340,700,800,300,400,3000,3100,440,815,1455,1000,2750,2770,2790,2810,2830,2850,1350,1400,1450,1500,700,800,1000,1100,1350,1400,900,1100,510.148,3050.56],'cm^-1')),
        HinderedRotor(inertia=(0.561197,'amu*angstrom^2'), symmetry=1, barrier=(12.903,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.870359,'amu*angstrom^2'), symmetry=1, barrier=(20.0113,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.00303129,'amu*angstrom^2'), symmetry=1, barrier=(20.0081,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.188911,'amu*angstrom^2'), symmetry=1, barrier=(4.34344,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(3.58667,'amu*angstrom^2'), symmetry=1, barrier=(82.4645,'kJ/mol'), semiclassical=False),
    ],
    spinMultiplicity = 2,
    opticalIsomers = 1,
    molecularWeight = (85.1674,'amu'),
    collisionModel = TransportData(shapeIndex=2, epsilon=(3225.18,'J/mol'), sigma=(6.08349,'angstroms'), dipoleMoment=(0,'C*m'), polarizability=(0,'angstroms^3'), rotrelaxcollnum=0, comment="""Epsilon & sigma estimated with Tc=503.77 K, Pc=32.5 bar (from Joback method)"""),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[0.896546,0.0589609,-1.64039e-05,-1.32731e-08,7.66956e-12,1024.55,26.7349], Tmin=(100,'K'), Tmax=(1044.35,'K')), NASAPolynomial(coeffs=[10.769,0.039543,-1.49348e-05,2.6552e-09,-1.80912e-13,-2040.65,-26.1253], Tmin=(1044.35,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(7.5235,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(436.51,'J/(mol*K)'), comment="""Thermo group additivity estimation: group(Cs-CsCsCsH) + group(Cs-CsCsHH) + longDistanceInteraction_noncyclic(CsCs-ST) + group(Cs-CsCsHH) + group(Cs-CsHHH) + group(Cs-CsHHH) + group(Cs-CsHHH) + radical(Isobutyl)"""),
)

species(
    label = 'CCC[C](C)C(92)',
    structure = SMILES('CCC[C](C)C'),
    E0 = (-12.1362,'kJ/mol'),
    modes = [
        HarmonicOscillator(frequencies=([2750,2762.5,2775,2787.5,2800,2812.5,2825,2837.5,2850,1350,1380,1410,1440,1470,1500,700,750,800,1000,1050,1100,1350,1375,1400,900,1000,1100,2750,2783.33,2816.67,2850,1425,1450,1225,1275,1270,1340,700,800,300,400,360,370,350,278.797,2890.67],'cm^-1')),
        HinderedRotor(inertia=(0.128796,'amu*angstrom^2'), symmetry=1, barrier=(7.10375,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.12879,'amu*angstrom^2'), symmetry=1, barrier=(7.10369,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.128792,'amu*angstrom^2'), symmetry=1, barrier=(7.10378,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.00216881,'amu*angstrom^2'), symmetry=1, barrier=(0.119627,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.12879,'amu*angstrom^2'), symmetry=1, barrier=(7.10378,'kJ/mol'), semiclassical=False),
    ],
    spinMultiplicity = 2,
    opticalIsomers = 1,
    molecularWeight = (85.1674,'amu'),
    collisionModel = TransportData(shapeIndex=2, epsilon=(3225.18,'J/mol'), sigma=(6.08349,'angstroms'), dipoleMoment=(0,'C*m'), polarizability=(0,'angstroms^3'), rotrelaxcollnum=0, comment="""Epsilon & sigma estimated with Tc=503.77 K, Pc=32.5 bar (from Joback method)"""),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[1.52175,0.0581738,-2.78879e-05,5.99589e-09,-4.95484e-13,-1374.47,23.1444], Tmin=(100,'K'), Tmax=(2807.76,'K')), NASAPolynomial(coeffs=[27.8591,0.0206525,-7.84249e-06,1.23629e-09,-7.16882e-14,-16164.1,-131.108], Tmin=(2807.76,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(-12.1362,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(436.51,'J/(mol*K)'), comment="""Thermo group additivity estimation: group(Cs-CsCsCsH) + group(Cs-CsCsHH) + longDistanceInteraction_noncyclic(CsCs-ST) + group(Cs-CsCsHH) + group(Cs-CsHHH) + group(Cs-CsHHH) + group(Cs-CsHHH) + radical(Tertalkyl)"""),
)

species(
    label = '[CH2]CCC(C)C(95)',
    structure = SMILES('[CH2]CCC(C)C'),
    E0 = (7.68753,'kJ/mol'),
    modes = [
        HarmonicOscillator(frequencies=([1380,1390,370,380,2900,435,2750,2783.33,2816.67,2850,1425,1450,1225,1275,1270,1340,700,800,300,400,3000,3100,440,815,1455,1000,2750,2770,2790,2810,2830,2850,1350,1400,1450,1500,700,800,1000,1100,1350,1400,900,1100,180,180],'cm^-1')),
        HinderedRotor(inertia=(0.0111309,'amu*angstrom^2'), symmetry=1, barrier=(20.6939,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.165292,'amu*angstrom^2'), symmetry=1, barrier=(3.8004,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.165263,'amu*angstrom^2'), symmetry=1, barrier=(3.79971,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.899986,'amu*angstrom^2'), symmetry=1, barrier=(20.6924,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.00182252,'amu*angstrom^2'), symmetry=1, barrier=(20.693,'kJ/mol'), semiclassical=False),
    ],
    spinMultiplicity = 2,
    opticalIsomers = 1,
    molecularWeight = (85.1674,'amu'),
    collisionModel = TransportData(shapeIndex=2, epsilon=(3225.18,'J/mol'), sigma=(6.08349,'angstroms'), dipoleMoment=(0,'C*m'), polarizability=(0,'angstroms^3'), rotrelaxcollnum=0, comment="""Epsilon & sigma estimated with Tc=503.77 K, Pc=32.5 bar (from Joback method)"""),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[0.900105,0.0581262,-1.44273e-05,-1.39324e-08,7.24906e-12,1044.61,26.0381], Tmin=(100,'K'), Tmax=(1100.21,'K')), NASAPolynomial(coeffs=[11.2261,0.0398903,-1.58862e-05,2.90083e-09,-2.00096e-13,-2395.98,-30.075], Tmin=(1100.21,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(7.68753,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(436.51,'J/(mol*K)'), comment="""Thermo group additivity estimation: group(Cs-CsCsCsH) + group(Cs-CsCsHH) + longDistanceInteraction_noncyclic(CsCs-ST) + group(Cs-CsCsHH) + group(Cs-CsHHH) + group(Cs-CsHHH) + group(Cs-CsHHH) + radical(RCCJ)"""),
)

species(
    label = 'CC[CH]C(C)C(93)',
    structure = SMILES('CC[CH]C(C)C'),
    E0 = (-3.01669,'kJ/mol'),
    modes = [
        HarmonicOscillator(frequencies=([1380,1390,370,380,2900,435,2750,2850,1437.5,1250,1305,750,350,3025,407.5,1350,352.5,2750,2762.5,2775,2787.5,2800,2812.5,2825,2837.5,2850,1350,1380,1410,1440,1470,1500,700,750,800,1000,1050,1100,1350,1375,1400,900,1000,1100,1144.93,4000],'cm^-1')),
        HinderedRotor(inertia=(0.173179,'amu*angstrom^2'), symmetry=1, barrier=(3.98172,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.17318,'amu*angstrom^2'), symmetry=1, barrier=(3.98175,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.173184,'amu*angstrom^2'), symmetry=1, barrier=(3.98185,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.902618,'amu*angstrom^2'), symmetry=1, barrier=(20.753,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.902618,'amu*angstrom^2'), symmetry=1, barrier=(20.753,'kJ/mol'), semiclassical=False),
    ],
    spinMultiplicity = 2,
    opticalIsomers = 1,
    molecularWeight = (85.1674,'amu'),
    collisionModel = TransportData(shapeIndex=2, epsilon=(3225.18,'J/mol'), sigma=(6.08349,'angstroms'), dipoleMoment=(0,'C*m'), polarizability=(0,'angstroms^3'), rotrelaxcollnum=0, comment="""Epsilon & sigma estimated with Tc=503.77 K, Pc=32.5 bar (from Joback method)"""),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[1.03399,0.0562814,-1.51385e-05,-9.19127e-09,4.69246e-12,-248.772,26.3353], Tmin=(100,'K'), Tmax=(1201.65,'K')), NASAPolynomial(coeffs=[10.3799,0.04146,-1.69702e-05,3.10547e-09,-2.12708e-13,-3670.92,-25.3636], Tmin=(1201.65,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(-3.01669,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(436.51,'J/(mol*K)'), comment="""Thermo group additivity estimation: group(Cs-CsCsCsH) + group(Cs-CsCsHH) + longDistanceInteraction_noncyclic(CsCs-ST) + group(Cs-CsCsHH) + group(Cs-CsHHH) + group(Cs-CsHHH) + group(Cs-CsHHH) + radical(Cs_S)"""),
)

species(
    label = '[CH2]C(CC)CC(83)',
    structure = SMILES('[CH2]C(CC)CC'),
    E0 = (10.8707,'kJ/mol'),
    modes = [
        HarmonicOscillator(frequencies=([1380,1390,370,380,2900,435,2750,2783.33,2816.67,2850,1425,1450,1225,1275,1270,1340,700,800,300,400,3000,3100,440,815,1455,1000,2750,2770,2790,2810,2830,2850,1350,1400,1450,1500,700,800,1000,1100,1350,1400,900,1100,510.148,3050.56],'cm^-1')),
        HinderedRotor(inertia=(0.561197,'amu*angstrom^2'), symmetry=1, barrier=(12.903,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.870359,'amu*angstrom^2'), symmetry=1, barrier=(20.0113,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.00303129,'amu*angstrom^2'), symmetry=1, barrier=(20.0081,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.188911,'amu*angstrom^2'), symmetry=1, barrier=(4.34344,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(3.58667,'amu*angstrom^2'), symmetry=1, barrier=(82.4645,'kJ/mol'), semiclassical=False),
    ],
    spinMultiplicity = 2,
    opticalIsomers = 1,
    molecularWeight = (85.1674,'amu'),
    collisionModel = TransportData(shapeIndex=2, epsilon=(3225.18,'J/mol'), sigma=(6.08349,'angstroms'), dipoleMoment=(0,'C*m'), polarizability=(0,'angstroms^3'), rotrelaxcollnum=0, comment="""Epsilon & sigma estimated with Tc=503.77 K, Pc=32.5 bar (from Joback method)"""),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[0.896546,0.0589609,-1.64039e-05,-1.32731e-08,7.66956e-12,1427.13,26.0417], Tmin=(100,'K'), Tmax=(1044.35,'K')), NASAPolynomial(coeffs=[10.769,0.039543,-1.49348e-05,2.6552e-09,-1.80912e-13,-1638.08,-26.8185], Tmin=(1044.35,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(10.8707,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(436.51,'J/(mol*K)'), comment="""Thermo group additivity estimation: group(Cs-CsCsCsH) + group(Cs-CsCsHH) + longDistanceInteraction_noncyclic(CsCs-ST) + group(Cs-CsCsHH) + longDistanceInteraction_noncyclic(CsCs-ST) + group(Cs-CsHHH) + group(Cs-CsHHH) + group(Cs-CsHHH) + radical(Isobutyl)"""),
)

species(
    label = 'C[CH]CC(C)C(94)',
    structure = SMILES('C[CH]CC(C)C'),
    E0 = (-3.11255,'kJ/mol'),
    modes = [
        HarmonicOscillator(frequencies=([1380,1390,370,380,2900,435,2750,2850,1437.5,1250,1305,750,350,3025,407.5,1350,352.5,2750,2762.5,2775,2787.5,2800,2812.5,2825,2837.5,2850,1350,1380,1410,1440,1470,1500,700,750,800,1000,1050,1100,1350,1375,1400,900,1000,1100,305.605,1420.04],'cm^-1')),
        HinderedRotor(inertia=(0.00180503,'amu*angstrom^2'), symmetry=1, barrier=(0.119627,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.0867051,'amu*angstrom^2'), symmetry=1, barrier=(5.74634,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.0867194,'amu*angstrom^2'), symmetry=1, barrier=(5.74638,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.0866957,'amu*angstrom^2'), symmetry=1, barrier=(5.74659,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.0867106,'amu*angstrom^2'), symmetry=1, barrier=(5.74653,'kJ/mol'), semiclassical=False),
    ],
    spinMultiplicity = 2,
    opticalIsomers = 1,
    molecularWeight = (85.1674,'amu'),
    collisionModel = TransportData(shapeIndex=2, epsilon=(3225.18,'J/mol'), sigma=(6.08349,'angstroms'), dipoleMoment=(0,'C*m'), polarizability=(0,'angstroms^3'), rotrelaxcollnum=0, comment="""Epsilon & sigma estimated with Tc=503.77 K, Pc=32.5 bar (from Joback method)"""),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[1.03026,0.0577495,-2.29435e-05,4.02646e-10,1.32482e-12,-261.329,26.5349], Tmin=(100,'K'), Tmax=(1358.89,'K')), NASAPolynomial(coeffs=[10.3462,0.0409868,-1.62066e-05,2.87014e-09,-1.91135e-13,-3777.39,-24.8874], Tmin=(1358.89,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(-3.11255,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(436.51,'J/(mol*K)'), comment="""Thermo group additivity estimation: group(Cs-CsCsCsH) + group(Cs-CsCsHH) + longDistanceInteraction_noncyclic(CsCs-ST) + group(Cs-CsCsHH) + group(Cs-CsHHH) + group(Cs-CsHHH) + group(Cs-CsHHH) + radical(RCCJC)"""),
)

species(
    label = 'CC[C](C)CC(134)',
    structure = SMILES('CC[C](C)CC'),
    E0 = (-8.78897,'kJ/mol'),
    modes = [
        HarmonicOscillator(frequencies=([2750,2762.5,2775,2787.5,2800,2812.5,2825,2837.5,2850,1350,1380,1410,1440,1470,1500,700,750,800,1000,1050,1100,1350,1375,1400,900,1000,1100,2750,2783.33,2816.67,2850,1425,1450,1225,1275,1270,1340,700,800,300,400,360,370,350,278.797,2890.67],'cm^-1')),
        HinderedRotor(inertia=(0.128796,'amu*angstrom^2'), symmetry=1, barrier=(7.10375,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.12879,'amu*angstrom^2'), symmetry=1, barrier=(7.10369,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.128792,'amu*angstrom^2'), symmetry=1, barrier=(7.10378,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.00216881,'amu*angstrom^2'), symmetry=1, barrier=(0.119627,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.12879,'amu*angstrom^2'), symmetry=1, barrier=(7.10378,'kJ/mol'), semiclassical=False),
    ],
    spinMultiplicity = 2,
    opticalIsomers = 1,
    molecularWeight = (85.1674,'amu'),
    collisionModel = TransportData(shapeIndex=2, epsilon=(3225.18,'J/mol'), sigma=(6.08349,'angstroms'), dipoleMoment=(0,'C*m'), polarizability=(0,'angstroms^3'), rotrelaxcollnum=0, comment="""Epsilon & sigma estimated with Tc=503.77 K, Pc=32.5 bar (from Joback method)"""),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[1.52175,0.0581738,-2.78879e-05,5.99589e-09,-4.95484e-13,-971.891,23.1444], Tmin=(100,'K'), Tmax=(2807.76,'K')), NASAPolynomial(coeffs=[27.8591,0.0206525,-7.84249e-06,1.23629e-09,-7.16882e-14,-15761.5,-131.108], Tmin=(2807.76,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(-8.78897,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(436.51,'J/(mol*K)'), comment="""Thermo group additivity estimation: group(Cs-CsCsCsH) + group(Cs-CsCsHH) + longDistanceInteraction_noncyclic(CsCs-ST) + group(Cs-CsCsHH) + longDistanceInteraction_noncyclic(CsCs-ST) + group(Cs-CsHHH) + group(Cs-CsHHH) + group(Cs-CsHHH) + radical(Tertalkyl)"""),
)

species(
    label = 'C[CH]C(C)CC(129)',
    structure = SMILES('C[CH]C(C)CC'),
    E0 = (0.330507,'kJ/mol'),
    modes = [
        HarmonicOscillator(frequencies=([1380,1390,370,380,2900,435,2750,2850,1437.5,1250,1305,750,350,3025,407.5,1350,352.5,2750,2762.5,2775,2787.5,2800,2812.5,2825,2837.5,2850,1350,1380,1410,1440,1470,1500,700,750,800,1000,1050,1100,1350,1375,1400,900,1000,1100,1144.93,4000],'cm^-1')),
        HinderedRotor(inertia=(0.173179,'amu*angstrom^2'), symmetry=1, barrier=(3.98172,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.17318,'amu*angstrom^2'), symmetry=1, barrier=(3.98175,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.173184,'amu*angstrom^2'), symmetry=1, barrier=(3.98185,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.902618,'amu*angstrom^2'), symmetry=1, barrier=(20.753,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.902618,'amu*angstrom^2'), symmetry=1, barrier=(20.753,'kJ/mol'), semiclassical=False),
    ],
    spinMultiplicity = 2,
    opticalIsomers = 1,
    molecularWeight = (85.1674,'amu'),
    collisionModel = TransportData(shapeIndex=2, epsilon=(3225.18,'J/mol'), sigma=(6.08349,'angstroms'), dipoleMoment=(0,'C*m'), polarizability=(0,'angstroms^3'), rotrelaxcollnum=0, comment="""Epsilon & sigma estimated with Tc=503.77 K, Pc=32.5 bar (from Joback method)"""),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[1.03399,0.0562814,-1.51385e-05,-9.19127e-09,4.69246e-12,153.803,27.0285], Tmin=(100,'K'), Tmax=(1201.65,'K')), NASAPolynomial(coeffs=[10.3799,0.04146,-1.69702e-05,3.10547e-09,-2.12708e-13,-3268.34,-24.6705], Tmin=(1201.65,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(0.330507,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(436.51,'J/(mol*K)'), comment="""Thermo group additivity estimation: group(Cs-CsCsCsH) + group(Cs-CsCsHH) + longDistanceInteraction_noncyclic(CsCs-ST) + group(Cs-CsCsHH) + longDistanceInteraction_noncyclic(CsCs-ST) + group(Cs-CsHHH) + group(Cs-CsHHH) + group(Cs-CsHHH) + radical(Cs_S)"""),
)

species(
    label = '[CH2]C(C)(C)CC(113)',
    structure = SMILES('[CH2]C(C)(C)CC'),
    E0 = (-1.27349,'kJ/mol'),
    modes = [
        HarmonicOscillator(frequencies=([350,400,1190,1240,450,2750,2850,1437.5,1250,1305,750,350,3000,3100,440,815,1455,1000,2750,2762.5,2775,2787.5,2800,2812.5,2825,2837.5,2850,1350,1380,1410,1440,1470,1500,700,750,800,1000,1050,1100,1350,1375,1400,900,1000,1100,180],'cm^-1')),
        HinderedRotor(inertia=(0.0883421,'amu*angstrom^2'), symmetry=1, barrier=(12.23,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.088269,'amu*angstrom^2'), symmetry=1, barrier=(12.2301,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.531927,'amu*angstrom^2'), symmetry=1, barrier=(12.23,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.0199362,'amu*angstrom^2'), symmetry=1, barrier=(2.76077,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.531917,'amu*angstrom^2'), symmetry=1, barrier=(12.2298,'kJ/mol'), semiclassical=False),
    ],
    spinMultiplicity = 2,
    opticalIsomers = 1,
    molecularWeight = (85.1674,'amu'),
    collisionModel = TransportData(shapeIndex=2, epsilon=(3254,'J/mol'), sigma=(6.12704,'angstroms'), dipoleMoment=(0,'C*m'), polarizability=(0,'angstroms^3'), rotrelaxcollnum=0, comment="""Epsilon & sigma estimated with Tc=508.27 K, Pc=32.1 bar (from Joback method)"""),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[0.766465,0.0592752,-1.08159e-05,-2.15448e-08,1.06062e-11,-26.6672,24.4434], Tmin=(100,'K'), Tmax=(1063.81,'K')), NASAPolynomial(coeffs=[12.9433,0.0379425,-1.52155e-05,2.8197e-09,-1.97375e-13,-4001.1,-41.5593], Tmin=(1063.81,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(-1.27349,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(436.51,'J/(mol*K)'), comment="""Thermo group additivity estimation: group(Cs-CsCsCsCs) + group(Cs-CsCsHH) + longDistanceInteraction_noncyclic(CsCs-SQ) + group(Cs-CsHHH) + group(Cs-CsHHH) + group(Cs-CsHHH) + group(Cs-CsHHH) + radical(Neopentyl)"""),
)

species(
    label = '[CH2]CC(C)CC(135)',
    structure = SMILES('[CH2]CC(C)CC'),
    E0 = (11.0347,'kJ/mol'),
    modes = [
        HarmonicOscillator(frequencies=([1380,1390,370,380,2900,435,2750,2783.33,2816.67,2850,1425,1450,1225,1275,1270,1340,700,800,300,400,3000,3100,440,815,1455,1000,2750,2770,2790,2810,2830,2850,1350,1400,1450,1500,700,800,1000,1100,1350,1400,900,1100,180,180],'cm^-1')),
        HinderedRotor(inertia=(0.0111309,'amu*angstrom^2'), symmetry=1, barrier=(20.6939,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.165292,'amu*angstrom^2'), symmetry=1, barrier=(3.8004,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.165263,'amu*angstrom^2'), symmetry=1, barrier=(3.79971,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.899986,'amu*angstrom^2'), symmetry=1, barrier=(20.6924,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.00182252,'amu*angstrom^2'), symmetry=1, barrier=(20.693,'kJ/mol'), semiclassical=False),
    ],
    spinMultiplicity = 2,
    opticalIsomers = 1,
    molecularWeight = (85.1674,'amu'),
    collisionModel = TransportData(shapeIndex=2, epsilon=(3225.18,'J/mol'), sigma=(6.08349,'angstroms'), dipoleMoment=(0,'C*m'), polarizability=(0,'angstroms^3'), rotrelaxcollnum=0, comment="""Epsilon & sigma estimated with Tc=503.77 K, Pc=32.5 bar (from Joback method)"""),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[0.900105,0.0581262,-1.44273e-05,-1.39324e-08,7.24906e-12,1447.19,26.7312], Tmin=(100,'K'), Tmax=(1100.21,'K')), NASAPolynomial(coeffs=[11.2261,0.0398903,-1.58862e-05,2.90083e-09,-2.00096e-13,-1993.41,-29.3818], Tmin=(1100.21,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(11.0347,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(436.51,'J/(mol*K)'), comment="""Thermo group additivity estimation: group(Cs-CsCsCsH) + group(Cs-CsCsHH) + longDistanceInteraction_noncyclic(CsCs-ST) + group(Cs-CsCsHH) + longDistanceInteraction_noncyclic(CsCs-ST) + group(Cs-CsHHH) + group(Cs-CsHHH) + group(Cs-CsHHH) + radical(RCCJ)"""),
)

species(
    label = '[CH2]C(C)C(C)C(130)',
    structure = SMILES('[CH2]C(C)C(C)C'),
    E0 = (1.8292,'kJ/mol'),
    modes = [
        HarmonicOscillator(frequencies=([1380,1383.33,1386.67,1390,370,373.333,376.667,380,2800,3000,430,440,3000,3100,440,815,1455,1000,2750,2762.5,2775,2787.5,2800,2812.5,2825,2837.5,2850,1350,1380,1410,1440,1470,1500,700,750,800,1000,1050,1100,1350,1375,1400,900,1000,1100,254.036],'cm^-1')),
        HinderedRotor(inertia=(0.00261279,'amu*angstrom^2'), symmetry=1, barrier=(0.119627,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.0950793,'amu*angstrom^2'), symmetry=1, barrier=(4.35445,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.0950188,'amu*angstrom^2'), symmetry=1, barrier=(4.35462,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.00261102,'amu*angstrom^2'), symmetry=1, barrier=(0.119627,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.459127,'amu*angstrom^2'), symmetry=1, barrier=(21.0387,'kJ/mol'), semiclassical=False),
    ],
    spinMultiplicity = 2,
    opticalIsomers = 1,
    molecularWeight = (85.1674,'amu'),
    collisionModel = TransportData(shapeIndex=2, epsilon=(3248.64,'J/mol'), sigma=(6.11673,'angstroms'), dipoleMoment=(0,'C*m'), polarizability=(0,'angstroms^3'), rotrelaxcollnum=0, comment="""Epsilon & sigma estimated with Tc=507.43 K, Pc=32.21 bar (from Joback method)"""),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[0.97594,0.055668,-3.71128e-06,-2.88953e-08,1.38443e-11,338.401,25.4757], Tmin=(100,'K'), Tmax=(996.475,'K')), NASAPolynomial(coeffs=[11.3813,0.0382086,-1.40223e-05,2.48449e-09,-1.70432e-13,-2942.25,-30.7434], Tmin=(996.475,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(1.8292,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(436.51,'J/(mol*K)'), comment="""Thermo group additivity estimation: group(Cs-CsCsCsH) + longDistanceInteraction_noncyclic(CsCs-TT) + group(Cs-CsCsCsH) + longDistanceInteraction_noncyclic(CsCs-TT) + group(Cs-CsHHH) + group(Cs-CsHHH) + group(Cs-CsHHH) + group(Cs-CsHHH) + radical(Isobutyl)"""),
)

species(
    label = 'C[CH]C(C)(C)C(160)',
    structure = SMILES('C[CH]C(C)(C)C'),
    E0 = (-11.7294,'kJ/mol'),
    modes = [
        HarmonicOscillator(frequencies=([350,400,1190,1240,450,3025,407.5,1350,352.5,2750,2759.09,2768.18,2777.27,2786.36,2795.45,2804.55,2813.64,2822.73,2831.82,2840.91,2850,1350,1371.43,1392.86,1414.29,1435.71,1457.14,1478.57,1500,700,733.333,766.667,800,1000,1033.33,1066.67,1100,1350,1366.67,1383.33,1400,900,966.667,1033.33,1100,229.126],'cm^-1')),
        HinderedRotor(inertia=(0.00321114,'amu*angstrom^2'), symmetry=1, barrier=(0.119627,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.00321107,'amu*angstrom^2'), symmetry=1, barrier=(0.119627,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.348885,'amu*angstrom^2'), symmetry=1, barrier=(12.9972,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.348894,'amu*angstrom^2'), symmetry=1, barrier=(12.9972,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.348875,'amu*angstrom^2'), symmetry=1, barrier=(12.9972,'kJ/mol'), semiclassical=False),
    ],
    spinMultiplicity = 2,
    opticalIsomers = 1,
    molecularWeight = (85.1674,'amu'),
    collisionModel = TransportData(shapeIndex=2, epsilon=(3254,'J/mol'), sigma=(6.12704,'angstroms'), dipoleMoment=(0,'C*m'), polarizability=(0,'angstroms^3'), rotrelaxcollnum=0, comment="""Epsilon & sigma estimated with Tc=508.27 K, Pc=32.1 bar (from Joback method)"""),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[0.956827,0.0550521,-2.76395e-06,-2.64205e-08,1.13667e-11,-1291.19,23.5068], Tmin=(100,'K'), Tmax=(1092.75,'K')), NASAPolynomial(coeffs=[11.9773,0.0399145,-1.658e-05,3.11425e-09,-2.18865e-13,-5204.46,-37.5233], Tmin=(1092.75,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(-11.7294,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(436.51,'J/(mol*K)'), comment="""Thermo group additivity estimation: group(Cs-CsCsCsCs) + group(Cs-CsCsHH) + longDistanceInteraction_noncyclic(CsCs-SQ) + group(Cs-CsHHH) + group(Cs-CsHHH) + group(Cs-CsHHH) + group(Cs-CsHHH) + radical(Cs_S)"""),
)

species(
    label = '[CH2]CC(C)(C)C(161)',
    structure = SMILES('[CH2]CC(C)(C)C'),
    E0 = (-1.02521,'kJ/mol'),
    modes = [
        HarmonicOscillator(frequencies=([350,400,1190,1240,450,2750,2850,1437.5,1250,1305,750,350,3000,3100,440,815,1455,1000,2750,2762.5,2775,2787.5,2800,2812.5,2825,2837.5,2850,1350,1380,1410,1440,1470,1500,700,750,800,1000,1050,1100,1350,1375,1400,900,1000,1100,180],'cm^-1')),
        HinderedRotor(inertia=(0.596145,'amu*angstrom^2'), symmetry=1, barrier=(13.7066,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.596215,'amu*angstrom^2'), symmetry=1, barrier=(13.7081,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.0416063,'amu*angstrom^2'), symmetry=1, barrier=(13.7091,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.0416578,'amu*angstrom^2'), symmetry=1, barrier=(13.7062,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.197854,'amu*angstrom^2'), symmetry=1, barrier=(4.54906,'kJ/mol'), semiclassical=False),
    ],
    spinMultiplicity = 2,
    opticalIsomers = 1,
    molecularWeight = (85.1674,'amu'),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[0.816257,0.0569162,-1.84528e-06,-3.18113e-08,1.43551e-11,2.58411,23.2376], Tmin=(100,'K'), Tmax=(1044.58,'K')), NASAPolynomial(coeffs=[13.3335,0.0375556,-1.50717e-05,2.81457e-09,-1.98684e-13,-4171.25,-45.158], Tmin=(1044.58,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(-1.02521,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(436.51,'J/(mol*K)'), comment="""Thermo group additivity estimation: group(Cs-CsCsCsCs) + group(Cs-CsCsHH) + longDistanceInteraction_noncyclic(CsCs-SQ) + group(Cs-CsHHH) + group(Cs-CsHHH) + group(Cs-CsHHH) + group(Cs-CsHHH) + radical(RCCJ)"""),
)

species(
    label = '[CH2]C[CH]C(44)',
    structure = SMILES('[CH2]C[CH]C'),
    E0 = (255.389,'kJ/mol'),
    modes = [
        HarmonicOscillator(frequencies=([2750,2800,2850,1350,1500,750,1050,1375,1000,2750,2850,1437.5,1250,1305,750,350,3000,3100,440,815,1455,1000,3025,407.5,1350,352.5,1553.58],'cm^-1')),
        HinderedRotor(inertia=(0.00260968,'amu*angstrom^2'), symmetry=1, barrier=(4.49712,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.195014,'amu*angstrom^2'), symmetry=1, barrier=(4.48376,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.194148,'amu*angstrom^2'), symmetry=1, barrier=(4.46384,'kJ/mol'), semiclassical=False),
    ],
    spinMultiplicity = 3,
    opticalIsomers = 1,
    molecularWeight = (56.1063,'amu'),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[2.59222,0.0286124,-6.18638e-06,-3.09349e-09,1.21714e-12,30768.7,19.1945], Tmin=(100,'K'), Tmax=(1492.36,'K')), NASAPolynomial(coeffs=[6.22727,0.02574,-1.02051e-05,1.7867e-09,-1.17175e-13,28918.6,-2.36123], Tmin=(1492.36,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(255.389,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(270.22,'J/(mol*K)'), comment="""Thermo group additivity estimation: group(Cs-CsCsHH) + group(Cs-CsCsHH) + group(Cs-CsHHH) + group(Cs-CsHHH) + radical(RCCJC) + radical(RCCJ)"""),
)

species(
    label = 'C[CH2](6)',
    structure = SMILES('C[CH2]'),
    E0 = (108.526,'kJ/mol'),
    modes = [
        HarmonicOscillator(frequencies=([2750,2800,2850,1350,1500,750,1050,1375,1000,484.904,1048.65,2314.79,2321.56,2325.77],'cm^-1')),
        HinderedRotor(inertia=(0.000785195,'amu*angstrom^2'), symmetry=1, barrier=(0.119627,'kJ/mol'), semiclassical=False),
    ],
    spinMultiplicity = 2,
    opticalIsomers = 1,
    molecularWeight = (29.0611,'amu'),
    collisionModel = TransportData(shapeIndex=2, epsilon=(2097.75,'J/mol'), sigma=(4.302,'angstroms'), dipoleMoment=(0,'C*m'), polarizability=(0,'angstroms^3'), rotrelaxcollnum=1.5, comment="""GRI-Mech"""),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[3.82181,-0.00343334,5.09247e-05,-6.20197e-08,2.37067e-11,13066,7.61651], Tmin=(100,'K'), Tmax=(900.32,'K')), NASAPolynomial(coeffs=[5.15627,0.00943113,-1.8194e-06,2.21182e-10,-1.43469e-14,12064.1,-2.9113], Tmin=(900.32,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(108.526,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(153.818,'J/(mol*K)'), comment="""Thermo group additivity estimation: group(Cs-CsHHH) + group(Cs-CsHHH) + radical(CCJ)"""),
)

species(
    label = '[CH2]CCC(3)',
    structure = SMILES('[CH2]CCC'),
    E0 = (60.9423,'kJ/mol'),
    modes = [
        HarmonicOscillator(frequencies=([2750,2800,2850,1350,1500,750,1050,1375,1000,2750,2783.33,2816.67,2850,1425,1450,1225,1275,1270,1340,700,800,300,400,3000,3100,440,815,1455,1000,180],'cm^-1')),
        HinderedRotor(inertia=(0.00258459,'amu*angstrom^2'), symmetry=1, barrier=(29.3456,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.605429,'amu*angstrom^2'), symmetry=1, barrier=(13.92,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.0101647,'amu*angstrom^2'), symmetry=1, barrier=(7.86263,'kJ/mol'), semiclassical=False),
    ],
    spinMultiplicity = 2,
    opticalIsomers = 1,
    molecularWeight = (57.1143,'amu'),
    collisionModel = TransportData(shapeIndex=2, epsilon=(2968.28,'J/mol'), sigma=(5.176,'angstroms'), dipoleMoment=(0,'C*m'), polarizability=(0,'angstroms^3'), rotrelaxcollnum=1.0, comment="""GRI-Mech"""),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[2.34849,0.0291322,9.83701e-06,-2.60513e-08,1.01154e-11,7395.26,17.3201], Tmin=(100,'K'), Tmax=(1048.1,'K')), NASAPolynomial(coeffs=[7.21568,0.0270038,-1.06548e-05,1.95473e-09,-1.35979e-13,5471.66,-10.6993], Tmin=(1048.1,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(60.9423,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(295.164,'J/(mol*K)'), comment="""Thermo group additivity estimation: group(Cs-CsCsHH) + group(Cs-CsCsHH) + group(Cs-CsHHH) + group(Cs-CsHHH) + radical(RCCJ)"""),
)

species(
    label = 'C=C(17)',
    structure = SMILES('C=C'),
    E0 = (42.1493,'kJ/mol'),
    modes = [
        HarmonicOscillator(frequencies=([2950,3000,3050,3100,1330,1430,900,1050,1000,1050,1600,1700],'cm^-1')),
    ],
    spinMultiplicity = 1,
    opticalIsomers = 1,
    molecularWeight = (28.0532,'amu'),
    collisionModel = TransportData(shapeIndex=2, epsilon=(2334.71,'J/mol'), sigma=(3.971,'angstroms'), dipoleMoment=(0,'C*m'), polarizability=(0,'angstroms^3'), rotrelaxcollnum=1.5, comment="""GRI-Mech"""),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[3.97978,-0.00757601,5.52989e-05,-6.36243e-08,2.31777e-11,5077.46,4.04611], Tmin=(100,'K'), Tmax=(940.437,'K')), NASAPolynomial(coeffs=[5.20289,0.00782461,-2.12694e-06,3.79716e-10,-2.94692e-14,3936.32,-6.62351], Tmin=(940.437,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(42.1493,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(133.032,'J/(mol*K)'), comment="""Thermo group additivity estimation: group(Cds-CdsHH) + group(Cds-CdsHH)"""),
)

species(
    label = 'H(8)',
    structure = SMILES('[H]'),
    E0 = (211.805,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
    molecularWeight = (1.00794,'amu'),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[2.5,-2.38914e-13,3.12709e-16,-1.33367e-19,1.7499e-23,25474.2,-0.444973], Tmin=(100,'K'), Tmax=(4383.16,'K')), NASAPolynomial(coeffs=[2.50003,-3.04997e-08,1.01101e-11,-1.48797e-15,8.20356e-20,25474.2,-0.445191], Tmin=(4383.16,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(211.805,'kJ/mol'), Cp0=(20.7862,'J/(mol*K)'), CpInf=(20.7862,'J/(mol*K)'), label="""H""", comment="""Thermo library: primaryThermoLibrary"""),
)

species(
    label = 'CC=CCCC(64)',
    structure = SMILES('CC=CCCC'),
    E0 = (-76.3276,'kJ/mol'),
    modes = [
        HarmonicOscillator(frequencies=([2750,2770,2790,2810,2830,2850,1350,1400,1450,1500,700,800,1000,1100,1350,1400,900,1100,2750,2783.33,2816.67,2850,1425,1450,1225,1275,1270,1340,700,800,300,400,2995,3025,975,1000,1300,1375,400,500,1630,1680,251.409,251.504],'cm^-1')),
        HinderedRotor(inertia=(0.00266637,'amu*angstrom^2'), symmetry=1, barrier=(0.119627,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.210301,'amu*angstrom^2'), symmetry=1, barrier=(9.42802,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.210056,'amu*angstrom^2'), symmetry=1, barrier=(9.42664,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.209988,'amu*angstrom^2'), symmetry=1, barrier=(9.4265,'kJ/mol'), semiclassical=False),
    ],
    spinMultiplicity = 1,
    opticalIsomers = 1,
    molecularWeight = (84.1595,'amu'),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[1.18465,0.0525691,-9.78365e-06,-1.51011e-08,7.02466e-12,-9070.98,24.2409], Tmin=(100,'K'), Tmax=(1123.91,'K')), NASAPolynomial(coeffs=[10.2411,0.03865,-1.56476e-05,2.87451e-09,-1.9854e-13,-12263.3,-25.6546], Tmin=(1123.91,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(-76.3276,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(415.724,'J/(mol*K)'), comment="""Thermo group additivity estimation: group(Cs-CsCsHH) + group(Cs-(Cds-Cds)CsHH) + group(Cs-CsHHH) + group(Cs-(Cds-Cds)HHH) + group(Cds-CdsCsH) + group(Cds-CdsCsH)"""),
)

species(
    label = 'C=CCCCC(48)',
    structure = SMILES('C=CCCCC'),
    E0 = (-64.0823,'kJ/mol'),
    modes = [
        HarmonicOscillator(frequencies=([2750,2800,2850,1350,1500,750,1050,1375,1000,2950,3100,1380,975,1025,1650,2750,2770,2790,2810,2830,2850,1425,1437.5,1450,1225,1250,1275,1270,1305,1340,700,750,800,300,350,400,3010,987.5,1337.5,450,1655,180,593.002,4000],'cm^-1')),
        HinderedRotor(inertia=(0.79719,'amu*angstrom^2'), symmetry=1, barrier=(18.329,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.797662,'amu*angstrom^2'), symmetry=1, barrier=(18.3398,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.796723,'amu*angstrom^2'), symmetry=1, barrier=(18.3182,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.154653,'amu*angstrom^2'), symmetry=1, barrier=(3.55577,'kJ/mol'), semiclassical=False),
    ],
    spinMultiplicity = 1,
    opticalIsomers = 1,
    molecularWeight = (84.1595,'amu'),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[1.14824,0.0522686,-5.10148e-06,-2.2328e-08,1.00605e-11,-7595.8,24.807], Tmin=(100,'K'), Tmax=(1073.47,'K')), NASAPolynomial(coeffs=[10.9452,0.0375731,-1.50431e-05,2.77295e-09,-1.93021e-13,-10955.8,-29.0055], Tmin=(1073.47,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(-64.0823,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(415.724,'J/(mol*K)'), comment="""Thermo group additivity estimation: group(Cs-CsCsHH) + group(Cs-CsCsHH) + group(Cs-(Cds-Cds)CsHH) + group(Cs-CsHHH) + group(Cds-CdsCsH) + group(Cds-CdsHH)"""),
)

species(
    label = '[CH2]CC(4)',
    structure = SMILES('[CH2]CC'),
    E0 = (84.7226,'kJ/mol'),
    modes = [
        HarmonicOscillator(frequencies=([2750,2800,2850,1350,1500,750,1050,1375,1000,2750,2850,1437.5,1250,1305,750,350,3000,3100,440,815,1455,1000],'cm^-1')),
        HinderedRotor(inertia=(0.0754772,'amu*angstrom^2'), symmetry=1, barrier=(1.73537,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.0764711,'amu*angstrom^2'), symmetry=1, barrier=(1.75822,'kJ/mol'), semiclassical=False),
    ],
    spinMultiplicity = 2,
    opticalIsomers = 1,
    molecularWeight = (43.0877,'amu'),
    collisionModel = TransportData(shapeIndex=2, epsilon=(2218.31,'J/mol'), sigma=(4.982,'angstroms'), dipoleMoment=(0,'C*m'), polarizability=(0,'angstroms^3'), rotrelaxcollnum=1.0, comment="""GRI-Mech"""),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[3.09191,0.0132173,2.75845e-05,-3.90845e-08,1.43312e-11,10228.4,12.4058], Tmin=(100,'K'), Tmax=(995.414,'K')), NASAPolynomial(coeffs=[5.69434,0.0196033,-7.42047e-06,1.35882e-09,-9.5621e-14,8875.84,-4.32905], Tmin=(995.414,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(84.7226,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(224.491,'J/(mol*K)'), comment="""Thermo group additivity estimation: group(Cs-CsCsHH) + group(Cs-CsHHH) + group(Cs-CsHHH) + radical(RCCJ)"""),
)

species(
    label = 'C=CC(36)',
    structure = SMILES('C=CC'),
    E0 = (6.12372,'kJ/mol'),
    modes = [
        HarmonicOscillator(frequencies=([2750,2800,2850,1350,1500,750,1050,1375,1000,3010,987.5,1337.5,450,1655,2950,3100,1380,975,1025,1650],'cm^-1')),
        HinderedRotor(inertia=(0.597443,'amu*angstrom^2'), symmetry=1, barrier=(13.7364,'kJ/mol'), semiclassical=False),
    ],
    spinMultiplicity = 1,
    opticalIsomers = 1,
    molecularWeight = (42.0797,'amu'),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[3.30972,0.00827549,3.37697e-05,-4.39283e-08,1.58761e-11,767.478,9.64367], Tmin=(100,'K'), Tmax=(988.018,'K')), NASAPolynomial(coeffs=[5.41223,0.0172863,-6.5134e-06,1.20318e-09,-8.55888e-14,-503.256,-4.80259], Tmin=(988.018,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(6.12372,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(203.705,'J/(mol*K)'), comment="""Thermo group additivity estimation: group(Cs-(Cds-Cds)HHH) + group(Cds-CdsCsH) + group(Cds-CdsHH)"""),
)

species(
    label = '[CH2][CH]C(38)',
    structure = SMILES('[CH2][CH]C'),
    E0 = (279.046,'kJ/mol'),
    modes = [
        HarmonicOscillator(frequencies=([2750,2800,2850,1350,1500,750,1050,1375,1000,3025,407.5,1350,352.5,3000,3100,440,815,1455,1000],'cm^-1')),
        HinderedRotor(inertia=(0.00418545,'amu*angstrom^2'), symmetry=1, barrier=(6.91847,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.00418539,'amu*angstrom^2'), symmetry=1, barrier=(6.91842,'kJ/mol'), semiclassical=False),
    ],
    spinMultiplicity = 3,
    opticalIsomers = 1,
    molecularWeight = (42.0797,'amu'),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[3.25506,0.0137285,1.00538e-05,-1.4379e-08,4.38758e-12,33590.4,14.1736], Tmin=(100,'K'), Tmax=(1201.85,'K')), NASAPolynomial(coeffs=[3.74307,0.0203097,-8.40109e-06,1.53861e-09,-1.05138e-13,32880.5,9.26403], Tmin=(1201.85,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(279.046,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(199.547,'J/(mol*K)'), comment="""Thermo group additivity estimation: group(Cs-CsCsHH) + group(Cs-CsHHH) + group(Cs-CsHHH) + radical(RCCJ) + radical(CCJC)"""),
)

species(
    label = 'C[CH]C[CH]CC(65)',
    structure = SMILES('C[CH]C[CH]CC'),
    E0 = (197.04,'kJ/mol'),
    modes = [
        HarmonicOscillator(frequencies=([2750,2770,2790,2810,2830,2850,1350,1400,1450,1500,700,800,1000,1100,1350,1400,900,1100,2750,2783.33,2816.67,2850,1425,1450,1225,1275,1270,1340,700,800,300,400,3000,3050,390,425,1340,1360,335,370,448.43,685.708,4000],'cm^-1')),
        HinderedRotor(inertia=(0.00320283,'amu*angstrom^2'), symmetry=1, barrier=(0.119627,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.00320283,'amu*angstrom^2'), symmetry=1, barrier=(0.119627,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.00320283,'amu*angstrom^2'), symmetry=1, barrier=(0.119627,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.00320283,'amu*angstrom^2'), symmetry=1, barrier=(0.119627,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.00320283,'amu*angstrom^2'), symmetry=1, barrier=(0.119627,'kJ/mol'), semiclassical=False),
    ],
    spinMultiplicity = 3,
    opticalIsomers = 1,
    molecularWeight = (84.1595,'amu'),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[2.33888,0.0488737,-1.82503e-05,1.90014e-09,1.05402e-13,23746.3,25.3213], Tmin=(100,'K'), Tmax=(2279.05,'K')), NASAPolynomial(coeffs=[21.3482,0.0263123,-1.05107e-05,1.71586e-09,-1.02518e-13,12276.4,-88.2013], Tmin=(2279.05,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(197.04,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(411.566,'J/(mol*K)'), comment="""Thermo group additivity estimation: group(Cs-CsCsHH) + group(Cs-CsCsHH) + group(Cs-CsCsHH) + group(Cs-CsCsHH) + group(Cs-CsHHH) + group(Cs-CsHHH) + radical(RCCJCC) + radical(RCCJC)"""),
)

species(
    label = '[CH2]CC[CH]C(30)',
    structure = SMILES('[CH2]CC[CH]C'),
    E0 = (231.608,'kJ/mol'),
    modes = [
        HarmonicOscillator(frequencies=([2750,2800,2850,1350,1500,750,1050,1375,1000,2750,2783.33,2816.67,2850,1425,1450,1225,1275,1270,1340,700,800,300,400,3000,3100,440,815,1455,1000,3025,407.5,1350,352.5,2179.49,2179.5],'cm^-1')),
        HinderedRotor(inertia=(0.00315792,'amu*angstrom^2'), symmetry=1, barrier=(10.6448,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.462981,'amu*angstrom^2'), symmetry=1, barrier=(10.6449,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(2.11859,'amu*angstrom^2'), symmetry=1, barrier=(48.7106,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.0235914,'amu*angstrom^2'), symmetry=1, barrier=(7.9597,'kJ/mol'), semiclassical=False),
    ],
    spinMultiplicity = 3,
    opticalIsomers = 1,
    molecularWeight = (70.1329,'amu'),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[1.95927,0.043164,-1.89531e-05,3.28957e-09,-1.25694e-13,27931,23.7177], Tmin=(100,'K'), Tmax=(1936.27,'K')), NASAPolynomial(coeffs=[12.6504,0.0264642,-1.01886e-05,1.70857e-09,-1.07055e-13,22781.1,-37.5325], Tmin=(1936.27,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(231.608,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(340.893,'J/(mol*K)'), comment="""Thermo group additivity estimation: group(Cs-CsCsHH) + group(Cs-CsCsHH) + group(Cs-CsCsHH) + group(Cs-CsHHH) + group(Cs-CsHHH) + radical(RCCJC) + radical(RCCJ)"""),
)

species(
    label = '[CH3](11)',
    structure = SMILES('[CH3]'),
    E0 = (135.382,'kJ/mol'),
    modes = [
        HarmonicOscillator(frequencies=([570.572,1407.44,1409.17,4000,4000,4000],'cm^-1')),
    ],
    spinMultiplicity = 2,
    opticalIsomers = 1,
    molecularWeight = (15.0345,'amu'),
    collisionModel = TransportData(shapeIndex=2, epsilon=(1197.29,'J/mol'), sigma=(3.8,'angstroms'), dipoleMoment=(0,'C*m'), polarizability=(0,'angstroms^3'), rotrelaxcollnum=0.0, comment="""GRI-Mech"""),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[3.91547,0.0018415,3.48754e-06,-3.32761e-09,8.50001e-13,16285.6,0.351726], Tmin=(100,'K'), Tmax=(1337.6,'K')), NASAPolynomial(coeffs=[3.5414,0.00476796,-1.82153e-06,3.28887e-10,-2.22554e-14,16224,1.6607], Tmin=(1337.6,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(135.382,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(108.088,'J/(mol*K)'), comment="""Thermo library: primaryThermoLibrary + radical(CH3)"""),
)

species(
    label = 'C[CH]CC[CH]C(66)',
    structure = SMILES('C[CH]CC[CH]C'),
    E0 = (197.028,'kJ/mol'),
    modes = [
        HarmonicOscillator(frequencies=([2750,2770,2790,2810,2830,2850,1350,1400,1450,1500,700,800,1000,1100,1350,1400,900,1100,2750,2783.33,2816.67,2850,1425,1450,1225,1275,1270,1340,700,800,300,400,3000,3050,390,425,1340,1360,335,370,200.056,728.032,3621.78],'cm^-1')),
        HinderedRotor(inertia=(0.00441921,'amu*angstrom^2'), symmetry=1, barrier=(0.119627,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.00441921,'amu*angstrom^2'), symmetry=1, barrier=(0.119627,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.00441921,'amu*angstrom^2'), symmetry=1, barrier=(0.119627,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.00441921,'amu*angstrom^2'), symmetry=1, barrier=(0.119627,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.00441921,'amu*angstrom^2'), symmetry=1, barrier=(0.119627,'kJ/mol'), semiclassical=False),
    ],
    spinMultiplicity = 3,
    opticalIsomers = 1,
    molecularWeight = (84.1595,'amu'),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[2.09537,0.0507271,-2.08289e-05,3.23712e-09,-1.13506e-13,23756.8,25.6821], Tmin=(100,'K'), Tmax=(2347.29,'K')), NASAPolynomial(coeffs=[21.8559,0.0251945,-9.71523e-06,1.55826e-09,-9.20701e-14,12237.3,-91.2897], Tmin=(2347.29,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(197.028,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(411.566,'J/(mol*K)'), comment="""Thermo group additivity estimation: group(Cs-CsCsHH) + group(Cs-CsCsHH) + group(Cs-CsCsHH) + group(Cs-CsCsHH) + group(Cs-CsHHH) + group(Cs-CsHHH) + radical(RCCJC) + radical(RCCJC)"""),
)

species(
    label = 'C[CH][CH]CCC(67)',
    structure = SMILES('C[CH][CH]CCC'),
    E0 = (197.04,'kJ/mol'),
    modes = [
        HarmonicOscillator(frequencies=([2750,2770,2790,2810,2830,2850,1350,1400,1450,1500,700,800,1000,1100,1350,1400,900,1100,2750,2783.33,2816.67,2850,1425,1450,1225,1275,1270,1340,700,800,300,400,3000,3050,390,425,1340,1360,335,370,448.43,685.709,4000],'cm^-1')),
        HinderedRotor(inertia=(0.00320283,'amu*angstrom^2'), symmetry=1, barrier=(0.119627,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.00320283,'amu*angstrom^2'), symmetry=1, barrier=(0.119627,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.00320283,'amu*angstrom^2'), symmetry=1, barrier=(0.119627,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.00320283,'amu*angstrom^2'), symmetry=1, barrier=(0.119627,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.00320283,'amu*angstrom^2'), symmetry=1, barrier=(0.119627,'kJ/mol'), semiclassical=False),
    ],
    spinMultiplicity = 3,
    opticalIsomers = 1,
    molecularWeight = (84.1595,'amu'),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[2.33888,0.0488737,-1.82503e-05,1.90014e-09,1.05402e-13,23746.3,25.3213], Tmin=(100,'K'), Tmax=(2279.05,'K')), NASAPolynomial(coeffs=[21.3482,0.0263123,-1.05107e-05,1.71586e-09,-1.02518e-13,12276.4,-88.2015], Tmin=(2279.05,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(197.04,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(411.566,'J/(mol*K)'), comment="""Thermo group additivity estimation: group(Cs-CsCsHH) + group(Cs-CsCsHH) + group(Cs-CsCsHH) + group(Cs-CsCsHH) + group(Cs-CsHHH) + group(Cs-CsHHH) + radical(RCCJC) + radical(RCCJCC)"""),
)

species(
    label = '[CH2]CCC[CH]C(53)',
    structure = SMILES('[CH2]CCC[CH]C'),
    E0 = (207.828,'kJ/mol'),
    modes = [
        HarmonicOscillator(frequencies=([2750,2800,2850,1350,1500,750,1050,1375,1000,2750,2770,2790,2810,2830,2850,1425,1437.5,1450,1225,1250,1275,1270,1305,1340,700,750,800,300,350,400,3000,3100,440,815,1455,1000,3025,407.5,1350,352.5,196.011,408.561,2850.35],'cm^-1')),
        HinderedRotor(inertia=(0.0727004,'amu*angstrom^2'), symmetry=1, barrier=(1.79832,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.0727004,'amu*angstrom^2'), symmetry=1, barrier=(1.79832,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.0727004,'amu*angstrom^2'), symmetry=1, barrier=(1.79832,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.0727004,'amu*angstrom^2'), symmetry=1, barrier=(1.79832,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.0727004,'amu*angstrom^2'), symmetry=1, barrier=(1.79832,'kJ/mol'), semiclassical=False),
    ],
    spinMultiplicity = 3,
    opticalIsomers = 1,
    molecularWeight = (84.1595,'amu'),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[1.49208,0.0562119,-2.8069e-05,6.53705e-09,-6.02902e-13,25084.2,27.6084], Tmin=(100,'K'), Tmax=(2373.07,'K')), NASAPolynomial(coeffs=[16.6366,0.0306851,-1.1934e-05,2.00437e-09,-1.25401e-13,17896.2,-58.5434], Tmin=(2373.07,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(207.828,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(411.566,'J/(mol*K)'), comment="""Thermo group additivity estimation: group(Cs-CsCsHH) + group(Cs-CsCsHH) + group(Cs-CsCsHH) + group(Cs-CsCsHH) + group(Cs-CsHHH) + group(Cs-CsHHH) + radical(RCCJ) + radical(RCCJC)"""),
)

species(
    label = '[CH2][CH]CCCC(54)',
    structure = SMILES('[CH2][CH]CCCC'),
    E0 = (207.828,'kJ/mol'),
    modes = [
        HarmonicOscillator(frequencies=([2750,2800,2850,1350,1500,750,1050,1375,1000,2750,2770,2790,2810,2830,2850,1425,1437.5,1450,1225,1250,1275,1270,1305,1340,700,750,800,300,350,400,3000,3100,440,815,1455,1000,3025,407.5,1350,352.5,196.011,408.562,2850.35],'cm^-1')),
        HinderedRotor(inertia=(0.0727004,'amu*angstrom^2'), symmetry=1, barrier=(1.79833,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.0727004,'amu*angstrom^2'), symmetry=1, barrier=(1.79833,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.0727004,'amu*angstrom^2'), symmetry=1, barrier=(1.79833,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.0727004,'amu*angstrom^2'), symmetry=1, barrier=(1.79833,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.0727004,'amu*angstrom^2'), symmetry=1, barrier=(1.79833,'kJ/mol'), semiclassical=False),
    ],
    spinMultiplicity = 3,
    opticalIsomers = 1,
    molecularWeight = (84.1595,'amu'),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[1.49208,0.0562119,-2.80689e-05,6.53704e-09,-6.02899e-13,25084.2,27.6084], Tmin=(100,'K'), Tmax=(2373.14,'K')), NASAPolynomial(coeffs=[16.6369,0.0306847,-1.19339e-05,2.00433e-09,-1.25398e-13,17896,-58.5452], Tmin=(2373.14,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(207.828,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(411.566,'J/(mol*K)'), comment="""Thermo group additivity estimation: group(Cs-CsCsHH) + group(Cs-CsCsHH) + group(Cs-CsCsHH) + group(Cs-CsCsHH) + group(Cs-CsHHH) + group(Cs-CsHHH) + radical(RCCJC) + radical(RCCJ)"""),
)

species(
    label = 'CH2(S)(14)',
    structure = SMILES('[CH2]'),
    E0 = (419.091,'kJ/mol'),
    modes = [
        HarmonicOscillator(frequencies=([1369.93,2896.01,2896.02],'cm^-1')),
    ],
    spinMultiplicity = 1,
    opticalIsomers = 1,
    molecularWeight = (14.0266,'amu'),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[4.10264,-0.00144068,5.45067e-06,-3.58e-09,7.56186e-13,50400.6,-0.411763], Tmin=(100,'K'), Tmax=(1442.37,'K')), NASAPolynomial(coeffs=[2.62649,0.00394761,-1.49923e-06,2.54537e-10,-1.62954e-14,50691.7,6.7837], Tmin=(1442.37,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(419.091,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(58.2013,'J/(mol*K)'), label="""CH2(S)""", comment="""Thermo library: primaryThermoLibrary"""),
)

species(
    label = 'C[CH]CCC(25)',
    structure = SMILES('C[CH]CCC'),
    E0 = (26.362,'kJ/mol'),
    modes = [
        HarmonicOscillator(frequencies=([2750,2770,2790,2810,2830,2850,1350,1400,1450,1500,700,800,1000,1100,1350,1400,900,1100,2750,2783.33,2816.67,2850,1425,1450,1225,1275,1270,1340,700,800,300,400,3025,407.5,1350,352.5,2618.67,2618.68],'cm^-1')),
        HinderedRotor(inertia=(0.0964395,'amu*angstrom^2'), symmetry=1, barrier=(8.14887,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.0964417,'amu*angstrom^2'), symmetry=1, barrier=(8.14886,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.0964456,'amu*angstrom^2'), symmetry=1, barrier=(8.14878,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.44131,'amu*angstrom^2'), symmetry=1, barrier=(37.289,'kJ/mol'), semiclassical=False),
    ],
    spinMultiplicity = 2,
    opticalIsomers = 1,
    molecularWeight = (71.1408,'amu'),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[1.77406,0.0441716,-1.45096e-05,-1.38062e-09,1.20145e-12,3254.59,22.595], Tmin=(100,'K'), Tmax=(1461.28,'K')), NASAPolynomial(coeffs=[8.55175,0.0345021,-1.37027e-05,2.4114e-09,-1.59038e-13,325.327,-15.9195], Tmin=(1461.28,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(26.362,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(365.837,'J/(mol*K)'), comment="""Thermo group additivity estimation: group(Cs-CsCsHH) + group(Cs-CsCsHH) + group(Cs-CsCsHH) + group(Cs-CsHHH) + group(Cs-CsHHH) + radical(RCCJC)"""),
)

species(
    label = '[CH]C(20)',
    structure = SMILES('[CH]C'),
    E0 = (351.472,'kJ/mol'),
    modes = [
        HarmonicOscillator(frequencies=([2750,2800,2850,1350,1500,750,1050,1375,1000,431.535,1804.51],'cm^-1')),
        HinderedRotor(inertia=(0.000906354,'amu*angstrom^2'), symmetry=1, barrier=(0.119627,'kJ/mol'), semiclassical=False),
    ],
    spinMultiplicity = 3,
    opticalIsomers = 1,
    molecularWeight = (28.0532,'amu'),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[3.73284,-0.000244647,3.59194e-05,-4.44283e-08,1.65885e-11,42287.5,7.07803], Tmin=(100,'K'), Tmax=(940.483,'K')), NASAPolynomial(coeffs=[5.42972,0.00816765,-2.42527e-06,4.22634e-10,-3.09411e-14,41277.1,-4.67909], Tmin=(940.483,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(351.472,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(128.874,'J/(mol*K)'), comment="""Thermo group additivity estimation: group(Cs-CsHHH) + group(Cs-CsHHH) + radical(CCJ2_triplet)"""),
)

species(
    label = '[CH]CCCC(33)',
    structure = SMILES('[CH]CCCC'),
    E0 = (280.131,'kJ/mol'),
    modes = [
        HarmonicOscillator(frequencies=([2750,2800,2850,1350,1500,750,1050,1375,1000,2750,2770,2790,2810,2830,2850,1425,1437.5,1450,1225,1250,1275,1270,1305,1340,700,750,800,300,350,400,200,800,1066.67,1333.33,1600],'cm^-1')),
        HinderedRotor(inertia=(0.156089,'amu*angstrom^2'), symmetry=1, barrier=(3.5888,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.156089,'amu*angstrom^2'), symmetry=1, barrier=(3.5888,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.156089,'amu*angstrom^2'), symmetry=1, barrier=(3.5888,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.156089,'amu*angstrom^2'), symmetry=1, barrier=(3.5888,'kJ/mol'), semiclassical=False),
    ],
    spinMultiplicity = 3,
    opticalIsomers = 1,
    molecularWeight = (70.1329,'amu'),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[1.49126,0.0476291,-1.77616e-05,-4.75466e-09,3.6807e-12,33788.6,21.8621], Tmin=(100,'K'), Tmax=(1127.23,'K')), NASAPolynomial(coeffs=[10.1002,0.0302049,-1.20403e-05,2.19072e-09,-1.50455e-13,31013.9,-24.4007], Tmin=(1127.23,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(280.131,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(340.893,'J/(mol*K)'), comment="""Thermo group additivity estimation: group(Cs-CsCsHH) + group(Cs-CsCsHH) + group(Cs-CsCsHH) + group(Cs-CsHHH) + group(Cs-CsHHH) + radical(CCJ2_triplet)"""),
)

species(
    label = 'C[C]CCCC(69)',
    structure = SMILES('C[C]CCCC'),
    E0 = (256.351,'kJ/mol'),
    modes = [
        HarmonicOscillator(frequencies=([2750,2770,2790,2810,2830,2850,1350,1400,1450,1500,700,800,1000,1100,1350,1400,900,1100,2750,2770,2790,2810,2830,2850,1425,1437.5,1450,1225,1250,1275,1270,1305,1340,700,750,800,300,350,400,200,800,1200,1600],'cm^-1')),
        HinderedRotor(inertia=(0.156089,'amu*angstrom^2'), symmetry=1, barrier=(3.5888,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.156089,'amu*angstrom^2'), symmetry=1, barrier=(3.5888,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.156089,'amu*angstrom^2'), symmetry=1, barrier=(3.5888,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.156089,'amu*angstrom^2'), symmetry=1, barrier=(3.5888,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.156089,'amu*angstrom^2'), symmetry=1, barrier=(3.5888,'kJ/mol'), semiclassical=False),
    ],
    spinMultiplicity = 3,
    opticalIsomers = 1,
    molecularWeight = (84.1595,'amu'),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[0.55961,0.065392,-4.03687e-05,1.25377e-08,-1.57983e-12,30964.4,26.3798], Tmin=(100,'K'), Tmax=(1818.37,'K')), NASAPolynomial(coeffs=[16.4429,0.0304519,-1.15457e-05,1.97026e-09,-1.2694e-13,25188.1,-59.7449], Tmin=(1818.37,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(256.351,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(411.566,'J/(mol*K)'), comment="""Thermo group additivity estimation: group(Cs-CsCsHH) + group(Cs-CsCsHH) + group(Cs-CsCsHH) + group(Cs-CsCsHH) + group(Cs-CsHHH) + group(Cs-CsHHH) + radical(CCJ2_triplet)"""),
)

species(
    label = '[CH2]C[CH2](27)',
    structure = SMILES('[CH2]C[CH2]'),
    E0 = (289.969,'kJ/mol'),
    modes = [
        HarmonicOscillator(frequencies=([2750,2850,1437.5,1250,1305,750,350,3000,3033.33,3066.67,3100,415,465,780,850,1435,1475,900,1100],'cm^-1')),
        HinderedRotor(inertia=(0.00735507,'amu*angstrom^2'), symmetry=1, barrier=(6.35285,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.00735344,'amu*angstrom^2'), symmetry=1, barrier=(6.3511,'kJ/mol'), semiclassical=False),
    ],
    spinMultiplicity = 3,
    opticalIsomers = 1,
    molecularWeight = (42.0797,'amu'),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[3.13018,0.013936,1.71981e-05,-2.69371e-08,9.94931e-12,34911,13.362], Tmin=(100,'K'), Tmax=(1011.03,'K')), NASAPolynomial(coeffs=[5.48786,0.01731,-6.6529e-06,1.21645e-09,-8.50336e-14,33785.1,-1.24883], Tmin=(1011.03,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(289.969,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(199.547,'J/(mol*K)'), comment="""Thermo group additivity estimation: group(Cs-CsCsHH) + group(Cs-CsHHH) + group(Cs-CsHHH) + radical(RCCJ) + radical(RCCJ)"""),
)

species(
    label = '[CH2]CC[CH2](29)',
    structure = SMILES('[CH2]CC[CH2]'),
    E0 = (266.189,'kJ/mol'),
    modes = [
        HarmonicOscillator(frequencies=([2750,2783.33,2816.67,2850,1425,1450,1225,1275,1270,1340,700,800,300,400,3000,3033.33,3066.67,3100,415,465,780,850,1435,1475,900,1100,1262.28],'cm^-1')),
        HinderedRotor(inertia=(0.0852535,'amu*angstrom^2'), symmetry=1, barrier=(1.96015,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.0855538,'amu*angstrom^2'), symmetry=1, barrier=(1.96705,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.0854654,'amu*angstrom^2'), symmetry=1, barrier=(1.96502,'kJ/mol'), semiclassical=False),
    ],
    spinMultiplicity = 3,
    opticalIsomers = 1,
    molecularWeight = (56.1063,'amu'),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[2.37803,0.0299377,-7.70878e-07,-1.37455e-08,5.72562e-12,32078.3,18.3087], Tmin=(100,'K'), Tmax=(1093.21,'K')), NASAPolynomial(coeffs=[7.15006,0.0244836,-9.76147e-06,1.78352e-09,-1.23053e-13,30317.5,-8.42031], Tmin=(1093.21,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(266.189,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(270.22,'J/(mol*K)'), comment="""Thermo group additivity estimation: group(Cs-CsCsHH) + group(Cs-CsCsHH) + group(Cs-CsHHH) + group(Cs-CsHHH) + radical(RCCJ) + radical(RCCJ)"""),
)

species(
    label = '[CH2]CC[CH]CC(51)',
    structure = SMILES('[CH2]CC[CH]CC'),
    E0 = (207.84,'kJ/mol'),
    modes = [
        HarmonicOscillator(frequencies=([2750,2800,2850,1350,1500,750,1050,1375,1000,2750,2770,2790,2810,2830,2850,1425,1437.5,1450,1225,1250,1275,1270,1305,1340,700,750,800,300,350,400,3000,3100,440,815,1455,1000,3025,407.5,1350,352.5,219.618,428.64,3500],'cm^-1')),
        HinderedRotor(inertia=(0.040347,'amu*angstrom^2'), symmetry=1, barrier=(1.30103,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.040347,'amu*angstrom^2'), symmetry=1, barrier=(1.30103,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.040347,'amu*angstrom^2'), symmetry=1, barrier=(1.30103,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.040347,'amu*angstrom^2'), symmetry=1, barrier=(1.30103,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.040347,'amu*angstrom^2'), symmetry=1, barrier=(1.30103,'kJ/mol'), semiclassical=False),
    ],
    spinMultiplicity = 3,
    opticalIsomers = 1,
    molecularWeight = (84.1595,'amu'),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[1.71179,0.0545483,-2.58693e-05,5.4579e-09,-4.38974e-13,25075.2,26.6482], Tmin=(100,'K'), Tmax=(2623.45,'K')), NASAPolynomial(coeffs=[25.9075,0.0202392,-7.72903e-06,1.22332e-09,-7.11968e-14,11491.4,-115.112], Tmin=(2623.45,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(207.84,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(411.566,'J/(mol*K)'), comment="""Thermo group additivity estimation: group(Cs-CsCsHH) + group(Cs-CsCsHH) + group(Cs-CsCsHH) + group(Cs-CsCsHH) + group(Cs-CsHHH) + group(Cs-CsHHH) + radical(RCCJCC) + radical(RCCJ)"""),
)

species(
    label = '[CH2][CH2](18)',
    structure = SMILES('[CH2][CH2]'),
    E0 = (313.796,'kJ/mol'),
    modes = [
        HarmonicOscillator(frequencies=([3000,3100,440,815,1455,1000,1603.88,1608.31,1610.47,1610.49,1612.55],'cm^-1')),
        HinderedRotor(inertia=(0.00570137,'amu*angstrom^2'), symmetry=1, barrier=(10.6016,'kJ/mol'), semiclassical=False),
    ],
    spinMultiplicity = 3,
    opticalIsomers = 1,
    molecularWeight = (28.0532,'amu'),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[3.86542,-0.0037011,4.71842e-05,-6.14997e-08,2.54414e-11,37752.3,8.63166], Tmin=(100,'K'), Tmax=(845.141,'K')), NASAPolynomial(coeffs=[5.89967,0.00441356,1.29138e-06,-4.58004e-10,3.6788e-14,36774.8,-4.58888], Tmin=(845.141,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(313.796,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(128.874,'J/(mol*K)'), comment="""Thermo group additivity estimation: group(Cs-CsHHH) + group(Cs-CsHHH) + radical(CCJ) + radical(CCJ)"""),
)

species(
    label = '[CH2]C[CH]CCC(52)',
    structure = SMILES('[CH2]C[CH]CCC'),
    E0 = (207.84,'kJ/mol'),
    modes = [
        HarmonicOscillator(frequencies=([2750,2800,2850,1350,1500,750,1050,1375,1000,2750,2770,2790,2810,2830,2850,1425,1437.5,1450,1225,1250,1275,1270,1305,1340,700,750,800,300,350,400,3000,3100,440,815,1455,1000,3025,407.5,1350,352.5,219.618,428.64,3500],'cm^-1')),
        HinderedRotor(inertia=(0.040347,'amu*angstrom^2'), symmetry=1, barrier=(1.30103,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.040347,'amu*angstrom^2'), symmetry=1, barrier=(1.30103,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.040347,'amu*angstrom^2'), symmetry=1, barrier=(1.30103,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.040347,'amu*angstrom^2'), symmetry=1, barrier=(1.30103,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.040347,'amu*angstrom^2'), symmetry=1, barrier=(1.30103,'kJ/mol'), semiclassical=False),
    ],
    spinMultiplicity = 3,
    opticalIsomers = 1,
    molecularWeight = (84.1595,'amu'),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[1.71179,0.0545483,-2.58693e-05,5.4579e-09,-4.38974e-13,25075.2,26.6482], Tmin=(100,'K'), Tmax=(2623.45,'K')), NASAPolynomial(coeffs=[25.9075,0.0202392,-7.72903e-06,1.22332e-09,-7.11968e-14,11491.4,-115.112], Tmin=(2623.45,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(207.84,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(411.566,'J/(mol*K)'), comment="""Thermo group additivity estimation: group(Cs-CsCsHH) + group(Cs-CsCsHH) + group(Cs-CsCsHH) + group(Cs-CsCsHH) + group(Cs-CsHHH) + group(Cs-CsHHH) + radical(RCCJCC) + radical(RCCJ)"""),
)

species(
    label = '[CH2]CCC[CH2](32)',
    structure = SMILES('[CH2]CCC[CH2]'),
    E0 = (242.408,'kJ/mol'),
    modes = [
        HarmonicOscillator(frequencies=([2750,2770,2790,2810,2830,2850,1425,1437.5,1450,1225,1250,1275,1270,1305,1340,700,750,800,300,350,400,3000,3033.33,3066.67,3100,415,465,780,850,1435,1475,900,1100,950.226,951.614],'cm^-1')),
        HinderedRotor(inertia=(0.151601,'amu*angstrom^2'), symmetry=1, barrier=(3.48561,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.152212,'amu*angstrom^2'), symmetry=1, barrier=(3.49966,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.151937,'amu*angstrom^2'), symmetry=1, barrier=(3.49333,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.152225,'amu*angstrom^2'), symmetry=1, barrier=(3.49995,'kJ/mol'), semiclassical=False),
    ],
    spinMultiplicity = 3,
    opticalIsomers = 1,
    molecularWeight = (70.1329,'amu'),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[1.59414,0.0462306,-1.93712e-05,-2.7787e-10,1.60866e-12,29247,23.3748], Tmin=(100,'K'), Tmax=(1249.03,'K')), NASAPolynomial(coeffs=[9.59758,0.0304327,-1.22077e-05,2.20147e-09,-1.4915e-13,26480.7,-20.0873], Tmin=(1249.03,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(242.408,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(340.893,'J/(mol*K)'), comment="""Thermo group additivity estimation: group(Cs-CsCsHH) + group(Cs-CsCsHH) + group(Cs-CsCsHH) + group(Cs-CsHHH) + group(Cs-CsHHH) + radical(RCCJ) + radical(RCCJ)"""),
)

species(
    label = '[CH2]CCCC[CH2](55)',
    structure = SMILES('[CH2]CCCC[CH2]'),
    E0 = (218.628,'kJ/mol'),
    modes = [
        HarmonicOscillator(frequencies=([2750,2764.29,2778.57,2792.86,2807.14,2821.43,2835.71,2850,1425,1433.33,1441.67,1450,1225,1241.67,1258.33,1275,1270,1293.33,1316.67,1340,700,733.333,766.667,800,300,333.333,366.667,400,3000,3033.33,3066.67,3100,415,465,780,850,1435,1475,900,1100,200,800,1600],'cm^-1')),
        HinderedRotor(inertia=(0.156089,'amu*angstrom^2'), symmetry=1, barrier=(3.5888,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.156089,'amu*angstrom^2'), symmetry=1, barrier=(3.5888,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.156089,'amu*angstrom^2'), symmetry=1, barrier=(3.5888,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.156089,'amu*angstrom^2'), symmetry=1, barrier=(3.5888,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.156089,'amu*angstrom^2'), symmetry=1, barrier=(3.5888,'kJ/mol'), semiclassical=False),
    ],
    spinMultiplicity = 3,
    opticalIsomers = 1,
    molecularWeight = (84.1595,'amu'),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[0.837596,0.0621036,-3.61211e-05,1.03934e-08,-1.21291e-12,26414.9,28.3503], Tmin=(100,'K'), Tmax=(1928.7,'K')), NASAPolynomial(coeffs=[15.7887,0.0310956,-1.20053e-05,2.05758e-09,-1.32395e-13,20647.7,-53.6008], Tmin=(1928.7,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(218.628,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(411.566,'J/(mol*K)'), comment="""Thermo group additivity estimation: group(Cs-CsCsHH) + group(Cs-CsCsHH) + group(Cs-CsCsHH) + group(Cs-CsCsHH) + group(Cs-CsHHH) + group(Cs-CsHHH) + radical(RCCJ) + radical(RCCJ)"""),
)

species(
    label = '[CH2]CCCC(5)',
    structure = SMILES('[CH2]CCCC'),
    E0 = (37.1621,'kJ/mol'),
    modes = [
        HarmonicOscillator(frequencies=([2750,2800,2850,1350,1500,750,1050,1375,1000,2750,2770,2790,2810,2830,2850,1425,1437.5,1450,1225,1250,1275,1270,1305,1340,700,750,800,300,350,400,3000,3100,440,815,1455,1000,183.76,4000],'cm^-1')),
        HinderedRotor(inertia=(0.004573,'amu*angstrom^2'), symmetry=1, barrier=(5.34159,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.222912,'amu*angstrom^2'), symmetry=1, barrier=(5.34157,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.922873,'amu*angstrom^2'), symmetry=1, barrier=(22.1142,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.922851,'amu*angstrom^2'), symmetry=1, barrier=(22.1142,'kJ/mol'), semiclassical=False),
    ],
    spinMultiplicity = 2,
    opticalIsomers = 1,
    molecularWeight = (71.1408,'amu'),
    collisionModel = TransportData(shapeIndex=2, epsilon=(3051.59,'J/mol'), sigma=(5.73385,'angstroms'), dipoleMoment=(0,'C*m'), polarizability=(0,'angstroms^3'), rotrelaxcollnum=0, comment="""Epsilon & sigma estimated with Tc=476.65 K, Pc=36.73 bar (from Joback method)"""),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[1.58666,0.0452266,-8.34997e-06,-1.27363e-08,5.90991e-12,4562.96,22.303], Tmin=(100,'K'), Tmax=(1123.46,'K')), NASAPolynomial(coeffs=[9.06519,0.033879,-1.35994e-05,2.48444e-09,-1.70983e-13,1918.36,-18.9386], Tmin=(1123.46,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(37.1621,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(365.837,'J/(mol*K)'), comment="""Thermo group additivity estimation: group(Cs-CsCsHH) + group(Cs-CsCsHH) + group(Cs-CsCsHH) + group(Cs-CsHHH) + group(Cs-CsHHH) + radical(RCCJ)"""),
)

species(
    label = 'CH2(T)(19)',
    structure = SMILES('[CH2]'),
    E0 = (381.37,'kJ/mol'),
    modes = [
        HarmonicOscillator(frequencies=([1066.91,2790.98,3622.38],'cm^-1')),
    ],
    spinMultiplicity = 3,
    opticalIsomers = 1,
    molecularWeight = (14.0266,'amu'),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[4.01192,-0.000154978,3.26298e-06,-2.40422e-09,5.69496e-13,45867.7,0.533201], Tmin=(100,'K'), Tmax=(1104.62,'K')), NASAPolynomial(coeffs=[3.14983,0.00296674,-9.76056e-07,1.54115e-10,-9.50339e-15,46058.1,4.77808], Tmin=(1104.62,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(381.37,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(58.2013,'J/(mol*K)'), label="""CH2(T)""", comment="""Thermo library: primaryThermoLibrary"""),
)

species(
    label = '[CH]CCCCC(56)',
    structure = SMILES('[CH]CCCCC'),
    E0 = (256.351,'kJ/mol'),
    modes = [
        HarmonicOscillator(frequencies=([2750,2800,2850,1350,1500,750,1050,1375,1000,2750,2764.29,2778.57,2792.86,2807.14,2821.43,2835.71,2850,1425,1433.33,1441.67,1450,1225,1241.67,1258.33,1275,1270,1293.33,1316.67,1340,700,733.333,766.667,800,300,333.333,366.667,400,200,800,1000,1200,1400,1600],'cm^-1')),
        HinderedRotor(inertia=(0.156089,'amu*angstrom^2'), symmetry=1, barrier=(3.5888,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.156089,'amu*angstrom^2'), symmetry=1, barrier=(3.5888,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.156089,'amu*angstrom^2'), symmetry=1, barrier=(3.5888,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.156089,'amu*angstrom^2'), symmetry=1, barrier=(3.5888,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.156089,'amu*angstrom^2'), symmetry=1, barrier=(3.5888,'kJ/mol'), semiclassical=False),
    ],
    spinMultiplicity = 3,
    opticalIsomers = 1,
    molecularWeight = (84.1595,'amu'),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[0.55961,0.065392,-4.03687e-05,1.25377e-08,-1.57983e-12,30964.4,27.4785], Tmin=(100,'K'), Tmax=(1818.37,'K')), NASAPolynomial(coeffs=[16.4429,0.0304519,-1.15457e-05,1.97026e-09,-1.2694e-13,25188.1,-58.6463], Tmin=(1818.37,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(256.351,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(411.566,'J/(mol*K)'), comment="""Thermo group additivity estimation: group(Cs-CsCsHH) + group(Cs-CsCsHH) + group(Cs-CsCsHH) + group(Cs-CsCsHH) + group(Cs-CsHHH) + group(Cs-CsHHH) + radical(CCJ2_triplet)"""),
)

species(
    label = 'CCC=CCC(81)',
    structure = SMILES('CCC=CCC'),
    E0 = (-75.1928,'kJ/mol'),
    modes = [
        HarmonicOscillator(frequencies=([2750,2770,2790,2810,2830,2850,1350,1400,1450,1500,700,800,1000,1100,1350,1400,900,1100,2750,2783.33,2816.67,2850,1425,1450,1225,1275,1270,1340,700,800,300,400,2995,3025,975,1000,1300,1375,400,500,1630,1680,298.811,299.372],'cm^-1')),
        HinderedRotor(inertia=(0.00187221,'amu*angstrom^2'), symmetry=1, barrier=(0.119627,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.162061,'amu*angstrom^2'), symmetry=1, barrier=(10.3332,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.162578,'amu*angstrom^2'), symmetry=1, barrier=(10.3534,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.160333,'amu*angstrom^2'), symmetry=1, barrier=(10.3684,'kJ/mol'), semiclassical=False),
    ],
    spinMultiplicity = 1,
    opticalIsomers = 1,
    molecularWeight = (84.1595,'amu'),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[1.27614,0.0486048,5.05475e-06,-3.27837e-08,1.37433e-11,-8935.87,23.8829], Tmin=(100,'K'), Tmax=(1047.22,'K')), NASAPolynomial(coeffs=[10.7488,0.0374726,-1.48806e-05,2.74931e-09,-1.92394e-13,-12293.4,-28.8123], Tmin=(1047.22,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(-75.1928,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(415.724,'J/(mol*K)'), comment="""Thermo group additivity estimation: group(Cs-(Cds-Cds)CsHH) + group(Cs-(Cds-Cds)CsHH) + group(Cs-CsHHH) + group(Cs-CsHHH) + group(Cds-CdsCsH) + group(Cds-CdsCsH)"""),
)

species(
    label = 'C=CCC(42)',
    structure = SMILES('C=CCC'),
    E0 = (-16.5218,'kJ/mol'),
    modes = [
        HarmonicOscillator(frequencies=([2750,2800,2850,1350,1500,750,1050,1375,1000,2950,3100,1380,975,1025,1650,2750,2850,1437.5,1250,1305,750,350,3010,987.5,1337.5,450,1655,385.428],'cm^-1')),
        HinderedRotor(inertia=(0.106004,'amu*angstrom^2'), symmetry=1, barrier=(11.0706,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.109803,'amu*angstrom^2'), symmetry=1, barrier=(11.1093,'kJ/mol'), semiclassical=False),
    ],
    spinMultiplicity = 1,
    opticalIsomers = 1,
    molecularWeight = (56.1063,'amu'),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[2.64174,0.0203702,3.05809e-05,-4.85557e-08,1.85231e-11,-1929.81,14.9537], Tmin=(100,'K'), Tmax=(991.795,'K')), NASAPolynomial(coeffs=[7.80836,0.0229241,-8.65886e-06,1.6005e-09,-1.13877e-13,-4105.11,-15.7295], Tmin=(991.795,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(-16.5218,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(274.378,'J/(mol*K)'), comment="""Thermo group additivity estimation: group(Cs-(Cds-Cds)CsHH) + group(Cs-CsHHH) + group(Cds-CdsCsH) + group(Cds-CdsHH)"""),
)

species(
    label = 'C=CCCC(24)',
    structure = SMILES('C=CCCC'),
    E0 = (-40.302,'kJ/mol'),
    modes = [
        HarmonicOscillator(frequencies=([2750,2800,2850,1350,1500,750,1050,1375,1000,2950,3100,1380,975,1025,1650,2750,2783.33,2816.67,2850,1425,1450,1225,1275,1270,1340,700,800,300,400,3010,987.5,1337.5,450,1655,433.779,435.457],'cm^-1')),
        HinderedRotor(inertia=(0.0746219,'amu*angstrom^2'), symmetry=1, barrier=(9.86906,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.0731864,'amu*angstrom^2'), symmetry=1, barrier=(9.85777,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.0733909,'amu*angstrom^2'), symmetry=1, barrier=(9.88638,'kJ/mol'), semiclassical=False),
    ],
    spinMultiplicity = 1,
    opticalIsomers = 1,
    molecularWeight = (70.1329,'amu'),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[1.90151,0.0362535,1.29138e-05,-3.55774e-08,1.43071e-11,-4763.1,19.8563], Tmin=(100,'K'), Tmax=(1027.61,'K')), NASAPolynomial(coeffs=[9.2806,0.0304043,-1.19376e-05,2.20665e-09,-1.55068e-13,-7487.39,-21.821], Tmin=(1027.61,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(-40.302,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(345.051,'J/(mol*K)'), comment="""Thermo group additivity estimation: group(Cs-CsCsHH) + group(Cs-(Cds-Cds)CsHH) + group(Cs-CsHHH) + group(Cds-CdsCsH) + group(Cds-CdsHH)"""),
)

species(
    label = '[CH2][CH]CC(45)',
    structure = SMILES('[CH2][CH]CC'),
    E0 = (255.389,'kJ/mol'),
    modes = [
        HarmonicOscillator(frequencies=([2750,2800,2850,1350,1500,750,1050,1375,1000,2750,2850,1437.5,1250,1305,750,350,3000,3100,440,815,1455,1000,3025,407.5,1350,352.5,1553.59],'cm^-1')),
        HinderedRotor(inertia=(0.00260965,'amu*angstrom^2'), symmetry=1, barrier=(4.49704,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.195017,'amu*angstrom^2'), symmetry=1, barrier=(4.48382,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.194149,'amu*angstrom^2'), symmetry=1, barrier=(4.46387,'kJ/mol'), semiclassical=False),
    ],
    spinMultiplicity = 3,
    opticalIsomers = 1,
    molecularWeight = (56.1063,'amu'),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[2.59222,0.0286124,-6.1863e-06,-3.09357e-09,1.21717e-12,30768.7,19.1945], Tmin=(100,'K'), Tmax=(1492.36,'K')), NASAPolynomial(coeffs=[6.2272,0.0257401,-1.02051e-05,1.78671e-09,-1.17176e-13,28918.6,-2.36081], Tmin=(1492.36,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(255.389,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(270.22,'J/(mol*K)'), comment="""Thermo group additivity estimation: group(Cs-CsCsHH) + group(Cs-CsCsHH) + group(Cs-CsHHH) + group(Cs-CsHHH) + radical(RCCJ) + radical(RCCJC)"""),
)

species(
    label = '[CH2]C[CH]CC(28)',
    structure = SMILES('[CH2]C[CH]CC'),
    E0 = (231.62,'kJ/mol'),
    modes = [
        HarmonicOscillator(frequencies=([2750,2800,2850,1350,1500,750,1050,1375,1000,2750,2783.33,2816.67,2850,1425,1450,1225,1275,1270,1340,700,800,300,400,3000,3100,440,815,1455,1000,3025,407.5,1350,352.5,1731.27,1731.45],'cm^-1')),
        HinderedRotor(inertia=(0.154565,'amu*angstrom^2'), symmetry=1, barrier=(3.55376,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.154503,'amu*angstrom^2'), symmetry=1, barrier=(3.55232,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.154543,'amu*angstrom^2'), symmetry=1, barrier=(3.55324,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.154529,'amu*angstrom^2'), symmetry=1, barrier=(3.55294,'kJ/mol'), semiclassical=False),
    ],
    spinMultiplicity = 3,
    opticalIsomers = 1,
    molecularWeight = (70.1329,'amu'),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[2.19946,0.0413401,-1.64432e-05,2.00795e-09,7.91308e-14,27920.7,22.6765], Tmin=(100,'K'), Tmax=(2000.48,'K')), NASAPolynomial(coeffs=[13.367,0.026096,-1.03257e-05,1.73981e-09,-1.08623e-13,22034.8,-42.4873], Tmin=(2000.48,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(231.62,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(340.893,'J/(mol*K)'), comment="""Thermo group additivity estimation: group(Cs-CsCsHH) + group(Cs-CsCsHH) + group(Cs-CsCsHH) + group(Cs-CsHHH) + group(Cs-CsHHH) + radical(RCCJ) + radical(RCCJCC)"""),
)

species(
    label = 'CC[CH][CH]CC(82)',
    structure = SMILES('CC[CH][CH]CC'),
    E0 = (197.052,'kJ/mol'),
    modes = [
        HarmonicOscillator(frequencies=([2750,2770,2790,2810,2830,2850,1350,1400,1450,1500,700,800,1000,1100,1350,1400,900,1100,2750,2783.33,2816.67,2850,1425,1450,1225,1275,1270,1340,700,800,300,400,3000,3050,390,425,1340,1360,335,370,414.175,1063.97,4000],'cm^-1')),
        HinderedRotor(inertia=(0.00298469,'amu*angstrom^2'), symmetry=1, barrier=(0.119627,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.00298469,'amu*angstrom^2'), symmetry=1, barrier=(0.119627,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.00298469,'amu*angstrom^2'), symmetry=1, barrier=(0.119627,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.00298469,'amu*angstrom^2'), symmetry=1, barrier=(0.119627,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.00298469,'amu*angstrom^2'), symmetry=1, barrier=(0.119627,'kJ/mol'), semiclassical=False),
    ],
    spinMultiplicity = 3,
    opticalIsomers = 1,
    molecularWeight = (84.1595,'amu'),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[2.5836,0.0470103,-1.56499e-05,5.46821e-10,3.2817e-13,23735.8,23.5696], Tmin=(100,'K'), Tmax=(2240.12,'K')), NASAPolynomial(coeffs=[21.1755,0.0270488,-1.11468e-05,1.84442e-09,-1.11017e-13,12085,-88.5337], Tmin=(2240.12,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(197.052,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(411.566,'J/(mol*K)'), comment="""Thermo group additivity estimation: group(Cs-CsCsHH) + group(Cs-CsCsHH) + group(Cs-CsCsHH) + group(Cs-CsCsHH) + group(Cs-CsHHH) + group(Cs-CsHHH) + radical(RCCJCC) + radical(RCCJCC)"""),
)

species(
    label = '[CH2][CH]CCC(31)',
    structure = SMILES('[CH2][CH]CCC'),
    E0 = (231.608,'kJ/mol'),
    modes = [
        HarmonicOscillator(frequencies=([2750,2800,2850,1350,1500,750,1050,1375,1000,2750,2783.33,2816.67,2850,1425,1450,1225,1275,1270,1340,700,800,300,400,3000,3100,440,815,1455,1000,3025,407.5,1350,352.5,2179.49,2179.5],'cm^-1')),
        HinderedRotor(inertia=(0.00315792,'amu*angstrom^2'), symmetry=1, barrier=(10.6448,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.462981,'amu*angstrom^2'), symmetry=1, barrier=(10.6449,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(2.11859,'amu*angstrom^2'), symmetry=1, barrier=(48.7106,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.0235914,'amu*angstrom^2'), symmetry=1, barrier=(7.9597,'kJ/mol'), semiclassical=False),
    ],
    spinMultiplicity = 3,
    opticalIsomers = 1,
    molecularWeight = (70.1329,'amu'),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[1.95927,0.043164,-1.89531e-05,3.28957e-09,-1.25694e-13,27931,23.7177], Tmin=(100,'K'), Tmax=(1936.27,'K')), NASAPolynomial(coeffs=[12.6504,0.0264642,-1.01886e-05,1.70857e-09,-1.07055e-13,22781.1,-37.5325], Tmin=(1936.27,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(231.608,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(340.893,'J/(mol*K)'), comment="""Thermo group additivity estimation: group(Cs-CsCsHH) + group(Cs-CsCsHH) + group(Cs-CsCsHH) + group(Cs-CsHHH) + group(Cs-CsHHH) + radical(RCCJC) + radical(RCCJ)"""),
)

species(
    label = 'CC[CH]CC(26)',
    structure = SMILES('CC[CH]CC'),
    E0 = (26.374,'kJ/mol'),
    modes = [
        HarmonicOscillator(frequencies=([2750,2770,2790,2810,2830,2850,1350,1400,1450,1500,700,800,1000,1100,1350,1400,900,1100,2750,2783.33,2816.67,2850,1425,1450,1225,1275,1270,1340,700,800,300,400,3025,407.5,1350,352.5,1749.16,1749.37],'cm^-1')),
        HinderedRotor(inertia=(0.10222,'amu*angstrom^2'), symmetry=1, barrier=(6.89316,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.102115,'amu*angstrom^2'), symmetry=1, barrier=(6.89875,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.102708,'amu*angstrom^2'), symmetry=1, barrier=(6.89428,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.102092,'amu*angstrom^2'), symmetry=1, barrier=(6.89724,'kJ/mol'), semiclassical=False),
    ],
    spinMultiplicity = 2,
    opticalIsomers = 1,
    molecularWeight = (71.1408,'amu'),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[1.9332,0.0431369,-1.41037e-05,-6.73514e-10,8.0467e-13,3248.32,21.1643], Tmin=(100,'K'), Tmax=(1648.17,'K')), NASAPolynomial(coeffs=[9.63726,0.0333951,-1.3388e-05,2.33364e-09,-1.51515e-13,-507.559,-23.5428], Tmin=(1648.17,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(26.374,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(365.837,'J/(mol*K)'), comment="""Thermo group additivity estimation: group(Cs-CsCsHH) + group(Cs-CsCsHH) + group(Cs-CsCsHH) + group(Cs-CsHHH) + group(Cs-CsHHH) + radical(RCCJCC)"""),
)

species(
    label = '[CH]CC(39)',
    structure = SMILES('[CH]CC'),
    E0 = (327.691,'kJ/mol'),
    modes = [
        HarmonicOscillator(frequencies=([2750,2800,2850,1350,1500,750,1050,1375,1000,2750,2850,1437.5,1250,1305,750,350,318.774,319.234,1779.44],'cm^-1')),
        HinderedRotor(inertia=(0.00165655,'amu*angstrom^2'), symmetry=1, barrier=(0.119627,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.00165208,'amu*angstrom^2'), symmetry=1, barrier=(0.119627,'kJ/mol'), semiclassical=False),
    ],
    spinMultiplicity = 3,
    opticalIsomers = 1,
    molecularWeight = (42.0797,'amu'),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[3.00517,0.0155111,1.85851e-05,-3.16862e-08,1.23732e-11,39453.6,11.9344], Tmin=(100,'K'), Tmax=(982.292,'K')), NASAPolynomial(coeffs=[6.73204,0.0159276,-5.86166e-06,1.06538e-09,-7.51285e-14,37969.2,-9.80821], Tmin=(982.292,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(327.691,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(199.547,'J/(mol*K)'), comment="""Thermo group additivity estimation: group(Cs-CsCsHH) + group(Cs-CsHHH) + group(Cs-CsHHH) + radical(CCJ2_triplet)"""),
)

species(
    label = '[CH]CCC(46)',
    structure = SMILES('[CH]CCC'),
    E0 = (303.911,'kJ/mol'),
    modes = [
        HarmonicOscillator(frequencies=([2750,2800,2850,1350,1500,750,1050,1375,1000,2750,2783.33,2816.67,2850,1425,1450,1225,1275,1270,1340,700,800,300,400,180,542.493,1300.89,2887.93],'cm^-1')),
        HinderedRotor(inertia=(0.0144352,'amu*angstrom^2'), symmetry=1, barrier=(17.3352,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.199942,'amu*angstrom^2'), symmetry=1, barrier=(4.59707,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.753961,'amu*angstrom^2'), symmetry=1, barrier=(17.3351,'kJ/mol'), semiclassical=False),
    ],
    spinMultiplicity = 3,
    opticalIsomers = 1,
    molecularWeight = (56.1063,'amu'),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[2.2609,0.0314424,7.48735e-07,-1.84873e-08,8.06281e-12,36620.5,16.8513], Tmin=(100,'K'), Tmax=(1038.41,'K')), NASAPolynomial(coeffs=[8.22231,0.0233776,-9.12316e-06,1.66749e-09,-1.15986e-13,34579.2,-16.0015], Tmin=(1038.41,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(303.911,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(270.22,'J/(mol*K)'), comment="""Thermo group additivity estimation: group(Cs-CsCsHH) + group(Cs-CsCsHH) + group(Cs-CsHHH) + group(Cs-CsHHH) + radical(CCJ2_triplet)"""),
)

species(
    label = 'CC[C]CCC(84)',
    structure = SMILES('CC[C]CCC'),
    E0 = (256.351,'kJ/mol'),
    modes = [
        HarmonicOscillator(frequencies=([2750,2770,2790,2810,2830,2850,1350,1400,1450,1500,700,800,1000,1100,1350,1400,900,1100,2750,2770,2790,2810,2830,2850,1425,1437.5,1450,1225,1250,1275,1270,1305,1340,700,750,800,300,350,400,200,800,1200,1600],'cm^-1')),
        HinderedRotor(inertia=(0.156089,'amu*angstrom^2'), symmetry=1, barrier=(3.5888,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.156089,'amu*angstrom^2'), symmetry=1, barrier=(3.5888,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.156089,'amu*angstrom^2'), symmetry=1, barrier=(3.5888,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.156089,'amu*angstrom^2'), symmetry=1, barrier=(3.5888,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.156089,'amu*angstrom^2'), symmetry=1, barrier=(3.5888,'kJ/mol'), semiclassical=False),
    ],
    spinMultiplicity = 3,
    opticalIsomers = 1,
    molecularWeight = (84.1595,'amu'),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[0.55961,0.065392,-4.03687e-05,1.25377e-08,-1.57983e-12,30964.4,26.3798], Tmin=(100,'K'), Tmax=(1818.37,'K')), NASAPolynomial(coeffs=[16.4429,0.0304519,-1.15457e-05,1.97026e-09,-1.2694e-13,25188.1,-59.7449], Tmin=(1818.37,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(256.351,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(411.566,'J/(mol*K)'), comment="""Thermo group additivity estimation: group(Cs-CsCsHH) + group(Cs-CsCsHH) + group(Cs-CsCsHH) + group(Cs-CsCsHH) + group(Cs-CsHHH) + group(Cs-CsHHH) + radical(CCJ2_triplet)"""),
)

species(
    label = 'C=C(C)CCC(91)',
    structure = SMILES('C=C(C)CCC'),
    E0 = (-79.357,'kJ/mol'),
    modes = [
        HarmonicOscillator(frequencies=([2750,2770,2790,2810,2830,2850,1350,1400,1450,1500,700,800,1000,1100,1350,1400,900,1100,2950,3100,1380,975,1025,1650,2750,2783.33,2816.67,2850,1425,1450,1225,1275,1270,1340,700,800,300,400,350,440,435,1725,316.942,316.945],'cm^-1')),
        HinderedRotor(inertia=(0.00167818,'amu*angstrom^2'), symmetry=1, barrier=(0.119627,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.139077,'amu*angstrom^2'), symmetry=1, barrier=(9.91399,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.139073,'amu*angstrom^2'), symmetry=1, barrier=(9.91398,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.139074,'amu*angstrom^2'), symmetry=1, barrier=(9.91399,'kJ/mol'), semiclassical=False),
    ],
    spinMultiplicity = 1,
    opticalIsomers = 1,
    molecularWeight = (84.1595,'amu'),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[0.999678,0.0568468,-1.92202e-05,-7.22608e-09,4.71876e-12,-9428.87,23.8557], Tmin=(100,'K'), Tmax=(1139.97,'K')), NASAPolynomial(coeffs=[11.0876,0.0375689,-1.50642e-05,2.74752e-09,-1.88774e-13,-12776.2,-30.7283], Tmin=(1139.97,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(-79.357,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(415.724,'J/(mol*K)'), comment="""Thermo group additivity estimation: group(Cs-CsCsHH) + group(Cs-(Cds-Cds)CsHH) + group(Cs-CsHHH) + group(Cs-(Cds-Cds)HHH) + group(Cds-CdsCsCs) + group(Cds-CdsHH)"""),
)

species(
    label = '[CH2][C](C)CCC(96)',
    structure = SMILES('[CH2][C](C)CCC'),
    E0 = (192.946,'kJ/mol'),
    modes = [
        HarmonicOscillator(frequencies=([2750,2770,2790,2810,2830,2850,1350,1400,1450,1500,700,800,1000,1100,1350,1400,900,1100,2750,2783.33,2816.67,2850,1425,1450,1225,1275,1270,1340,700,800,300,400,3000,3100,440,815,1455,1000,360,370,350,180,2476.06],'cm^-1')),
        HinderedRotor(inertia=(0.144223,'amu*angstrom^2'), symmetry=1, barrier=(3.31598,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.144225,'amu*angstrom^2'), symmetry=1, barrier=(3.31601,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.144238,'amu*angstrom^2'), symmetry=1, barrier=(3.31631,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.144204,'amu*angstrom^2'), symmetry=1, barrier=(3.31554,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.144245,'amu*angstrom^2'), symmetry=1, barrier=(3.31648,'kJ/mol'), semiclassical=False),
    ],
    spinMultiplicity = 3,
    opticalIsomers = 1,
    molecularWeight = (84.1595,'amu'),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[2.03401,0.0509563,1.07262e-05,-9.44335e-08,8.02874e-11,23269,23.9518], Tmin=(100,'K'), Tmax=(471.198,'K')), NASAPolynomial(coeffs=[3.66017,0.0490308,-2.09591e-05,3.89767e-09,-2.69542e-13,22983.9,15.9309], Tmin=(471.198,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(192.946,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(411.566,'J/(mol*K)'), comment="""Thermo group additivity estimation: group(Cs-CsCsCsH) + group(Cs-CsCsHH) + longDistanceInteraction_noncyclic(CsCs-ST) + group(Cs-CsCsHH) + group(Cs-CsHHH) + group(Cs-CsHHH) + group(Cs-CsHHH) + radical(Tertalkyl) + radical(Isobutyl)"""),
)

species(
    label = '[CH2]C([CH2])C(97)',
    structure = SMILES('[CH2]C([CH2])C'),
    E0 = (256.819,'kJ/mol'),
    modes = [
        HarmonicOscillator(frequencies=([1380,1390,370,380,2900,435,3000,3033.33,3066.67,3100,415,465,780,850,1435,1475,900,1100,2750,2800,2850,1350,1500,750,1050,1375,1000],'cm^-1')),
        HinderedRotor(inertia=(0.127339,'amu*angstrom^2'), symmetry=1, barrier=(2.92777,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.127458,'amu*angstrom^2'), symmetry=1, barrier=(2.9305,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.0224565,'amu*angstrom^2'), symmetry=1, barrier=(69.4846,'kJ/mol'), semiclassical=False),
    ],
    spinMultiplicity = 3,
    opticalIsomers = 1,
    molecularWeight = (56.1063,'amu'),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[2.45366,0.028127,9.38804e-06,-3.11347e-08,1.46677e-11,30949.4,17.7469], Tmin=(100,'K'), Tmax=(890.509,'K')), NASAPolynomial(coeffs=[7.5673,0.021302,-6.31007e-06,9.76125e-10,-6.24427e-14,29398.5,-9.92464], Tmin=(890.509,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(256.819,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(270.22,'J/(mol*K)'), comment="""Thermo group additivity estimation: group(Cs-CsCsCsH) + group(Cs-CsHHH) + group(Cs-CsHHH) + group(Cs-CsHHH) + radical(Isobutyl) + radical(Isobutyl)"""),
)

species(
    label = '[CH2]C(C)[CH]CC(98)',
    structure = SMILES('[CH2]C(C)[CH]CC'),
    E0 = (202.066,'kJ/mol'),
    modes = [
        HarmonicOscillator(frequencies=([1380,1390,370,380,2900,435,3000,3100,440,815,1455,1000,2750,2850,1437.5,1250,1305,750,350,3025,407.5,1350,352.5,2750,2770,2790,2810,2830,2850,1350,1400,1450,1500,700,800,1000,1100,1350,1400,900,1100,180,1249.07],'cm^-1')),
        HinderedRotor(inertia=(0.0904743,'amu*angstrom^2'), symmetry=1, barrier=(2.08018,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.0912646,'amu*angstrom^2'), symmetry=1, barrier=(2.09835,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.0900749,'amu*angstrom^2'), symmetry=1, barrier=(2.071,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.0908681,'amu*angstrom^2'), symmetry=1, barrier=(2.08924,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.0904328,'amu*angstrom^2'), symmetry=1, barrier=(2.07923,'kJ/mol'), semiclassical=False),
    ],
    spinMultiplicity = 3,
    opticalIsomers = 1,
    molecularWeight = (84.1595,'amu'),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[0.867061,0.0598704,-3.31987e-05,9.27409e-09,-1.05724e-12,24423.4,29.4302], Tmin=(100,'K'), Tmax=(1969.66,'K')), NASAPolynomial(coeffs=[15.0165,0.0311353,-1.13151e-05,1.8671e-09,-1.17089e-13,18849.6,-48.4237], Tmin=(1969.66,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(202.066,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(411.566,'J/(mol*K)'), comment="""Thermo group additivity estimation: group(Cs-CsCsCsH) + group(Cs-CsCsHH) + longDistanceInteraction_noncyclic(CsCs-ST) + group(Cs-CsCsHH) + group(Cs-CsHHH) + group(Cs-CsHHH) + group(Cs-CsHHH) + radical(Isobutyl) + radical(Cs_S)"""),
)

species(
    label = '[CH2]CC([CH2])C(99)',
    structure = SMILES('[CH2]CC([CH2])C'),
    E0 = (236.55,'kJ/mol'),
    modes = [
        HarmonicOscillator(frequencies=([1380,1390,370,380,2900,435,2750,2850,1437.5,1250,1305,750,350,3000,3033.33,3066.67,3100,415,465,780,850,1435,1475,900,1100,2750,2800,2850,1350,1500,750,1050,1375,1000,2661.19],'cm^-1')),
        HinderedRotor(inertia=(0.63139,'amu*angstrom^2'), symmetry=1, barrier=(14.5169,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.13674,'amu*angstrom^2'), symmetry=1, barrier=(3.14393,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.00361991,'amu*angstrom^2'), symmetry=1, barrier=(3.2004,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(3.42024,'amu*angstrom^2'), symmetry=1, barrier=(78.6381,'kJ/mol'), semiclassical=False),
    ],
    spinMultiplicity = 3,
    opticalIsomers = 1,
    molecularWeight = (70.1329,'amu'),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[1.69267,0.0435924,-8.4332e-06,-1.49552e-08,7.84384e-12,28539.7,23.4182], Tmin=(100,'K'), Tmax=(1005.54,'K')), NASAPolynomial(coeffs=[8.98465,0.0299451,-1.09882e-05,1.92993e-09,-1.31045e-13,26296.7,-15.6629], Tmin=(1005.54,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(236.55,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(340.893,'J/(mol*K)'), comment="""Thermo group additivity estimation: group(Cs-CsCsCsH) + group(Cs-CsCsHH) + longDistanceInteraction_noncyclic(CsCs-ST) + group(Cs-CsHHH) + group(Cs-CsHHH) + group(Cs-CsHHH) + radical(Isobutyl) + radical(RCCJ)"""),
)

species(
    label = '[CH2]C(C)C[CH]C(100)',
    structure = SMILES('[CH2]C(C)C[CH]C'),
    E0 = (201.97,'kJ/mol'),
    modes = [
        HarmonicOscillator(frequencies=([1380,1390,370,380,2900,435,3000,3100,440,815,1455,1000,2750,2850,1437.5,1250,1305,750,350,3025,407.5,1350,352.5,2750,2770,2790,2810,2830,2850,1350,1400,1450,1500,700,800,1000,1100,1350,1400,900,1100,180,1489.26],'cm^-1')),
        HinderedRotor(inertia=(0.0785987,'amu*angstrom^2'), symmetry=1, barrier=(1.80714,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.077129,'amu*angstrom^2'), symmetry=1, barrier=(1.77335,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.0780873,'amu*angstrom^2'), symmetry=1, barrier=(1.79538,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.0777716,'amu*angstrom^2'), symmetry=1, barrier=(1.78812,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.0773859,'amu*angstrom^2'), symmetry=1, barrier=(1.77926,'kJ/mol'), semiclassical=False),
    ],
    spinMultiplicity = 3,
    opticalIsomers = 1,
    molecularWeight = (84.1595,'amu'),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[1.19756,0.0579365,-3.12855e-05,8.83564e-09,-1.05773e-12,24394.8,28.3897], Tmin=(100,'K'), Tmax=(1786.76,'K')), NASAPolynomial(coeffs=[10.4272,0.0372741,-1.39393e-05,2.36346e-09,-1.5215e-13,21096.6,-21.4948], Tmin=(1786.76,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(201.97,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(411.566,'J/(mol*K)'), comment="""Thermo group additivity estimation: group(Cs-CsCsCsH) + group(Cs-CsCsHH) + longDistanceInteraction_noncyclic(CsCs-ST) + group(Cs-CsCsHH) + group(Cs-CsHHH) + group(Cs-CsHHH) + group(Cs-CsHHH) + radical(RCCJC) + radical(Isobutyl)"""),
)

species(
    label = '[CH2]C([CH2])CCC(101)',
    structure = SMILES('[CH2]C([CH2])CCC'),
    E0 = (212.606,'kJ/mol'),
    modes = [
        HarmonicOscillator(frequencies=([1380,1390,370,380,2900,435,2750,2783.33,2816.67,2850,1425,1450,1225,1275,1270,1340,700,800,300,400,3000,3033.33,3066.67,3100,415,465,780,850,1435,1475,900,1100,2750,2800,2850,1350,1500,750,1050,1375,1000,262.787,2072.95],'cm^-1')),
        HinderedRotor(inertia=(0.00236384,'amu*angstrom^2'), symmetry=1, barrier=(0.119627,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.0952473,'amu*angstrom^2'), symmetry=1, barrier=(5.06336,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.00232619,'amu*angstrom^2'), symmetry=1, barrier=(0.119783,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.0946509,'amu*angstrom^2'), symmetry=1, barrier=(5.03324,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(1.41248,'amu*angstrom^2'), symmetry=1, barrier=(70.6583,'kJ/mol'), semiclassical=False),
    ],
    spinMultiplicity = 3,
    opticalIsomers = 1,
    molecularWeight = (84.1595,'amu'),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[0.94629,0.0603019,-2.78544e-05,-1.90136e-09,4.43183e-12,25686.5,27.643], Tmin=(100,'K'), Tmax=(992.809,'K')), NASAPolynomial(coeffs=[10.2432,0.036692,-1.31042e-05,2.24241e-09,-1.49184e-13,23158.1,-20.5791], Tmin=(992.809,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(212.606,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(411.566,'J/(mol*K)'), comment="""Thermo group additivity estimation: group(Cs-CsCsCsH) + group(Cs-CsCsHH) + longDistanceInteraction_noncyclic(CsCs-ST) + group(Cs-CsCsHH) + group(Cs-CsHHH) + group(Cs-CsHHH) + group(Cs-CsHHH) + radical(Isobutyl) + radical(Isobutyl)"""),
)

species(
    label = '[CH2]CCC([CH2])C(102)',
    structure = SMILES('[CH2]CCC([CH2])C'),
    E0 = (212.77,'kJ/mol'),
    modes = [
        HarmonicOscillator(frequencies=([1380,1390,370,380,2900,435,2750,2783.33,2816.67,2850,1425,1450,1225,1275,1270,1340,700,800,300,400,3000,3033.33,3066.67,3100,415,465,780,850,1435,1475,900,1100,2750,2800,2850,1350,1500,750,1050,1375,1000,1075.22,2799.79],'cm^-1')),
        HinderedRotor(inertia=(0.10897,'amu*angstrom^2'), symmetry=1, barrier=(2.50544,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.107115,'amu*angstrom^2'), symmetry=1, barrier=(2.46279,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.691222,'amu*angstrom^2'), symmetry=1, barrier=(15.8926,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.108477,'amu*angstrom^2'), symmetry=1, barrier=(2.4941,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(3.51947,'amu*angstrom^2'), symmetry=1, barrier=(80.9194,'kJ/mol'), semiclassical=False),
    ],
    spinMultiplicity = 3,
    opticalIsomers = 1,
    molecularWeight = (84.1595,'amu'),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[0.905828,0.0600018,-2.7815e-05,3.76906e-11,2.86971e-12,25708.5,28.4894], Tmin=(100,'K'), Tmax=(1094.78,'K')), NASAPolynomial(coeffs=[10.8091,0.0368474,-1.39422e-05,2.46082e-09,-1.66083e-13,22759.3,-23.7512], Tmin=(1094.78,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(212.77,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(411.566,'J/(mol*K)'), comment="""Thermo group additivity estimation: group(Cs-CsCsCsH) + group(Cs-CsCsHH) + longDistanceInteraction_noncyclic(CsCs-ST) + group(Cs-CsCsHH) + group(Cs-CsHHH) + group(Cs-CsHHH) + group(Cs-CsHHH) + radical(RCCJ) + radical(Isobutyl)"""),
)

species(
    label = '[CH2]C(C)CC(103)',
    structure = SMILES('[CH2]C(C)CC'),
    E0 = (31.3037,'kJ/mol'),
    modes = [
        HarmonicOscillator(frequencies=([1380,1390,370,380,2900,435,2750,2850,1437.5,1250,1305,750,350,3000,3100,440,815,1455,1000,2750,2770,2790,2810,2830,2850,1350,1400,1450,1500,700,800,1000,1100,1350,1400,900,1100,3368.84],'cm^-1')),
        HinderedRotor(inertia=(0.0202244,'amu*angstrom^2'), symmetry=1, barrier=(13.7994,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.210977,'amu*angstrom^2'), symmetry=1, barrier=(4.85077,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(3.58286,'amu*angstrom^2'), symmetry=1, barrier=(82.3771,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.599675,'amu*angstrom^2'), symmetry=1, barrier=(13.7877,'kJ/mol'), semiclassical=False),
    ],
    spinMultiplicity = 2,
    opticalIsomers = 1,
    molecularWeight = (71.1408,'amu'),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[1.65934,0.0428139,2.16912e-06,-2.7391e-08,1.23521e-11,3856.87,21.7513], Tmin=(100,'K'), Tmax=(991.726,'K')), NASAPolynomial(coeffs=[9.17437,0.0322672,-1.17724e-05,2.07626e-09,-1.41961e-13,1394.37,-19.3421], Tmin=(991.726,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(31.3037,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(365.837,'J/(mol*K)'), comment="""Thermo group additivity estimation: group(Cs-CsCsCsH) + group(Cs-CsCsHH) + longDistanceInteraction_noncyclic(CsCs-ST) + group(Cs-CsHHH) + group(Cs-CsHHH) + group(Cs-CsHHH) + radical(Isobutyl)"""),
)

species(
    label = '[CH]C(C)CCC(104)',
    structure = SMILES('[CH]C(C)CCC'),
    E0 = (250.656,'kJ/mol'),
    modes = [
        HarmonicOscillator(frequencies=([1380,1390,370,380,2900,435,2750,2783.33,2816.67,2850,1425,1450,1225,1275,1270,1340,700,800,300,400,2750,2770,2790,2810,2830,2850,1350,1400,1450,1500,700,800,1000,1100,1350,1400,900,1100,200,800,1066.67,1333.33,1600],'cm^-1')),
        HinderedRotor(inertia=(0.156089,'amu*angstrom^2'), symmetry=1, barrier=(3.5888,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.156089,'amu*angstrom^2'), symmetry=1, barrier=(3.5888,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.156089,'amu*angstrom^2'), symmetry=1, barrier=(3.5888,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.156089,'amu*angstrom^2'), symmetry=1, barrier=(3.5888,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.156089,'amu*angstrom^2'), symmetry=1, barrier=(3.5888,'kJ/mol'), semiclassical=False),
    ],
    spinMultiplicity = 3,
    opticalIsomers = 1,
    molecularWeight = (84.1595,'amu'),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[0.808226,0.0604888,-2.37075e-05,-6.10858e-09,5.08143e-12,30270.1,26.2777], Tmin=(100,'K'), Tmax=(1098.6,'K')), NASAPolynomial(coeffs=[12.2373,0.0362549,-1.43486e-05,2.6121e-09,-1.79974e-13,26710.1,-34.7092], Tmin=(1098.6,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(250.656,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(411.566,'J/(mol*K)'), comment="""Thermo group additivity estimation: group(Cs-CsCsCsH) + group(Cs-CsCsHH) + longDistanceInteraction_noncyclic(CsCs-ST) + group(Cs-CsCsHH) + group(Cs-CsHHH) + group(Cs-CsHHH) + group(Cs-CsHHH) + radical(CCJ2_triplet)"""),
)

species(
    label = 'CCC=C(C)C(105)',
    structure = SMILES('CCC=C(C)C'),
    E0 = (-91.6024,'kJ/mol'),
    modes = [
        HarmonicOscillator(frequencies=([2750,2762.5,2775,2787.5,2800,2812.5,2825,2837.5,2850,1350,1380,1410,1440,1470,1500,700,750,800,1000,1050,1100,1350,1375,1400,900,1000,1100,3010,987.5,1337.5,450,1655,2750,2850,1437.5,1250,1305,750,350,350,440,435,1725,262.84],'cm^-1')),
        HinderedRotor(inertia=(0.174652,'amu*angstrom^2'), symmetry=1, barrier=(8.48771,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.177628,'amu*angstrom^2'), symmetry=1, barrier=(8.51283,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.174547,'amu*angstrom^2'), symmetry=1, barrier=(8.50099,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.177543,'amu*angstrom^2'), symmetry=1, barrier=(8.50725,'kJ/mol'), semiclassical=False),
    ],
    spinMultiplicity = 1,
    opticalIsomers = 1,
    molecularWeight = (84.1595,'amu'),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[1.02076,0.0572653,-2.40349e-05,-1.96769e-10,1.90923e-12,-10903.3,23.3489], Tmin=(100,'K'), Tmax=(1260.33,'K')), NASAPolynomial(coeffs=[11.0038,0.0376957,-1.51618e-05,2.73618e-09,-1.85344e-13,-14381.8,-30.9404], Tmin=(1260.33,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(-91.6024,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(415.724,'J/(mol*K)'), comment="""Thermo group additivity estimation: group(Cs-(Cds-Cds)CsHH) + group(Cs-CsHHH) + group(Cs-(Cds-Cds)HHH) + group(Cs-(Cds-Cds)HHH) + group(Cds-CdsCsCs) + group(Cds-CdsCsH)"""),
)

species(
    label = 'C=C(C)C(106)',
    structure = SMILES('C=C(C)C'),
    E0 = (-32.9313,'kJ/mol'),
    modes = [
        HarmonicOscillator(frequencies=([2750,2770,2790,2810,2830,2850,1350,1400,1450,1500,700,800,1000,1100,1350,1400,900,1100,2950,3100,1380,975,1025,1650,350,440,435,1725],'cm^-1')),
        HinderedRotor(inertia=(0.452593,'amu*angstrom^2'), symmetry=1, barrier=(10.406,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.453645,'amu*angstrom^2'), symmetry=1, barrier=(10.4302,'kJ/mol'), semiclassical=False),
    ],
    spinMultiplicity = 1,
    opticalIsomers = 1,
    molecularWeight = (56.1063,'amu'),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[2.44003,0.028538,2.55258e-06,-1.63976e-08,6.47075e-12,-3899.73,12.8316], Tmin=(100,'K'), Tmax=(1089.91,'K')), NASAPolynomial(coeffs=[6.76118,0.0251857,-1.00461e-05,1.83699e-09,-1.26819e-13,-5584.49,-11.7952], Tmin=(1089.91,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(-32.9313,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(274.378,'J/(mol*K)'), comment="""Thermo group additivity estimation: group(Cs-(Cds-Cds)HHH) + group(Cs-(Cds-Cds)HHH) + group(Cds-CdsCsCs) + group(Cds-CdsHH)"""),
)

species(
    label = '[CH2][C](C)C(107)',
    structure = SMILES('[CH2][C](C)C'),
    E0 = (237.159,'kJ/mol'),
    modes = [
        HarmonicOscillator(frequencies=([2750,2770,2790,2810,2830,2850,1350,1400,1450,1500,700,800,1000,1100,1350,1400,900,1100,3000,3100,440,815,1455,1000,360,370,350],'cm^-1')),
        HinderedRotor(inertia=(0.00324468,'amu*angstrom^2'), symmetry=1, barrier=(8.01922,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.00322008,'amu*angstrom^2'), symmetry=1, barrier=(7.97596,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.347292,'amu*angstrom^2'), symmetry=1, barrier=(7.98492,'kJ/mol'), semiclassical=False),
    ],
    spinMultiplicity = 3,
    opticalIsomers = 1,
    molecularWeight = (56.1063,'amu'),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[2.93003,0.0294162,-1.08085e-05,1.44449e-09,-2.49174e-14,28556.2,15.3611], Tmin=(100,'K'), Tmax=(2439.2,'K')), NASAPolynomial(coeffs=[14.0561,0.0160592,-5.60072e-06,8.42776e-10,-4.74577e-14,21674.2,-51.2177], Tmin=(2439.2,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(237.159,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(270.22,'J/(mol*K)'), comment="""Thermo group additivity estimation: group(Cs-CsCsCsH) + group(Cs-CsHHH) + group(Cs-CsHHH) + group(Cs-CsHHH) + radical(Tertalkyl) + radical(Isobutyl)"""),
)

species(
    label = '[CH2]C[C](C)C(108)',
    structure = SMILES('[CH2]C[C](C)C'),
    E0 = (216.89,'kJ/mol'),
    modes = [
        HarmonicOscillator(frequencies=([2750,2770,2790,2810,2830,2850,1350,1400,1450,1500,700,800,1000,1100,1350,1400,900,1100,2750,2850,1437.5,1250,1305,750,350,3000,3100,440,815,1455,1000,360,370,350,2672.94],'cm^-1')),
        HinderedRotor(inertia=(0.152554,'amu*angstrom^2'), symmetry=1, barrier=(3.50752,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.152576,'amu*angstrom^2'), symmetry=1, barrier=(3.50802,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.152526,'amu*angstrom^2'), symmetry=1, barrier=(3.50688,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.152573,'amu*angstrom^2'), symmetry=1, barrier=(3.50796,'kJ/mol'), semiclassical=False),
    ],
    spinMultiplicity = 3,
    opticalIsomers = 1,
    molecularWeight = (70.1329,'amu'),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[2.38036,0.0422226,-1.85277e-05,3.31325e-09,-1.89675e-13,26137.6,19.5923], Tmin=(100,'K'), Tmax=(2401.03,'K')), NASAPolynomial(coeffs=[19.9683,0.0185871,-7.30118e-06,1.1788e-09,-6.97534e-14,16058.8,-84.0653], Tmin=(2401.03,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(216.89,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(340.893,'J/(mol*K)'), comment="""Thermo group additivity estimation: group(Cs-CsCsCsH) + group(Cs-CsCsHH) + longDistanceInteraction_noncyclic(CsCs-ST) + group(Cs-CsHHH) + group(Cs-CsHHH) + group(Cs-CsHHH) + radical(RCCJ) + radical(Tertalkyl)"""),
)

species(
    label = 'C[CH]C[C](C)C(109)',
    structure = SMILES('C[CH]C[C](C)C'),
    E0 = (182.31,'kJ/mol'),
    modes = [
        HarmonicOscillator(frequencies=([2750,2762.5,2775,2787.5,2800,2812.5,2825,2837.5,2850,1350,1380,1410,1440,1470,1500,700,750,800,1000,1050,1100,1350,1375,1400,900,1000,1100,2750,2850,1437.5,1250,1305,750,350,3025,407.5,1350,352.5,360,370,350,2265.11,2268.61],'cm^-1')),
        HinderedRotor(inertia=(0.19478,'amu*angstrom^2'), symmetry=1, barrier=(4.47838,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.194763,'amu*angstrom^2'), symmetry=1, barrier=(4.47797,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.194849,'amu*angstrom^2'), symmetry=1, barrier=(4.47996,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.196375,'amu*angstrom^2'), symmetry=1, barrier=(4.51504,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.194907,'amu*angstrom^2'), symmetry=1, barrier=(4.48128,'kJ/mol'), semiclassical=False),
    ],
    spinMultiplicity = 3,
    opticalIsomers = 1,
    molecularWeight = (84.1595,'amu'),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[1.47197,0.0613147,-5.81854e-05,5.05575e-08,-2.03514e-11,22012.3,26.0648], Tmin=(100,'K'), Tmax=(791.088,'K')), NASAPolynomial(coeffs=[0.703326,0.0543099,-2.4252e-05,4.5577e-09,-3.14691e-13,22474.7,31.7469], Tmin=(791.088,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(182.31,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(411.566,'J/(mol*K)'), comment="""Thermo group additivity estimation: group(Cs-CsCsCsH) + group(Cs-CsCsHH) + longDistanceInteraction_noncyclic(CsCs-ST) + group(Cs-CsCsHH) + group(Cs-CsHHH) + group(Cs-CsHHH) + group(Cs-CsHHH) + radical(Tertalkyl) + radical(RCCJC)"""),
)

species(
    label = 'CC[CH][C](C)C(110)',
    structure = SMILES('CC[CH][C](C)C'),
    E0 = (182.406,'kJ/mol'),
    modes = [
        HarmonicOscillator(frequencies=([2750,2762.5,2775,2787.5,2800,2812.5,2825,2837.5,2850,1350,1380,1410,1440,1470,1500,700,750,800,1000,1050,1100,1350,1375,1400,900,1000,1100,2750,2850,1437.5,1250,1305,750,350,3025,407.5,1350,352.5,360,370,350,180,3774.43],'cm^-1')),
        HinderedRotor(inertia=(0.0961369,'amu*angstrom^2'), symmetry=1, barrier=(2.21038,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.09613,'amu*angstrom^2'), symmetry=1, barrier=(2.21022,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.096127,'amu*angstrom^2'), symmetry=1, barrier=(2.21015,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.0961443,'amu*angstrom^2'), symmetry=1, barrier=(2.21055,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.0961466,'amu*angstrom^2'), symmetry=1, barrier=(2.2106,'kJ/mol'), semiclassical=False),
    ],
    spinMultiplicity = 3,
    opticalIsomers = 1,
    molecularWeight = (84.1595,'amu'),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[2.19765,0.0516106,-2.20518e-05,3.48251e-09,-1.22645e-13,21991.5,23.2461], Tmin=(100,'K'), Tmax=(2362.29,'K')), NASAPolynomial(coeffs=[25.0431,0.0221863,-9.24759e-06,1.52827e-09,-9.1428e-14,8614.47,-112.077], Tmin=(2362.29,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(182.406,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(411.566,'J/(mol*K)'), comment="""Thermo group additivity estimation: group(Cs-CsCsCsH) + group(Cs-CsCsHH) + longDistanceInteraction_noncyclic(CsCs-ST) + group(Cs-CsCsHH) + group(Cs-CsHHH) + group(Cs-CsHHH) + group(Cs-CsHHH) + radical(Cs_S) + radical(Tertalkyl)"""),
)

species(
    label = '[CH2]CC[C](C)C(111)',
    structure = SMILES('[CH2]CC[C](C)C'),
    E0 = (193.11,'kJ/mol'),
    modes = [
        HarmonicOscillator(frequencies=([2750,2770,2790,2810,2830,2850,1350,1400,1450,1500,700,800,1000,1100,1350,1400,900,1100,2750,2783.33,2816.67,2850,1425,1450,1225,1275,1270,1340,700,800,300,400,3000,3100,440,815,1455,1000,360,370,350,180,2549.89],'cm^-1')),
        HinderedRotor(inertia=(0.133245,'amu*angstrom^2'), symmetry=1, barrier=(3.06356,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.133245,'amu*angstrom^2'), symmetry=1, barrier=(3.06357,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.13324,'amu*angstrom^2'), symmetry=1, barrier=(3.06344,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.13324,'amu*angstrom^2'), symmetry=1, barrier=(3.06345,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.133244,'amu*angstrom^2'), symmetry=1, barrier=(3.06353,'kJ/mol'), semiclassical=False),
    ],
    spinMultiplicity = 3,
    opticalIsomers = 1,
    molecularWeight = (84.1595,'amu'),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[1.97979,0.0547295,-2.65263e-05,5.76402e-09,-4.86714e-13,23286.6,23.2217], Tmin=(100,'K'), Tmax=(2565.56,'K')), NASAPolynomial(coeffs=[17.1226,0.03112,-1.27225e-05,2.17706e-09,-1.37183e-13,15516.7,-64.101], Tmin=(2565.56,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(193.11,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(411.566,'J/(mol*K)'), comment="""Thermo group additivity estimation: group(Cs-CsCsCsH) + group(Cs-CsCsHH) + longDistanceInteraction_noncyclic(CsCs-ST) + group(Cs-CsCsHH) + group(Cs-CsHHH) + group(Cs-CsHHH) + group(Cs-CsHHH) + radical(Tertalkyl) + radical(RCCJ)"""),
)

species(
    label = 'CC[C](C)C(112)',
    structure = SMILES('CC[C](C)C'),
    E0 = (11.6441,'kJ/mol'),
    modes = [
        HarmonicOscillator(frequencies=([2750,2762.5,2775,2787.5,2800,2812.5,2825,2837.5,2850,1350,1380,1410,1440,1470,1500,700,750,800,1000,1050,1100,1350,1375,1400,900,1000,1100,2750,2850,1437.5,1250,1305,750,350,360,370,350,237.385],'cm^-1')),
        HinderedRotor(inertia=(0.00118783,'amu*angstrom^2'), symmetry=1, barrier=(6.31216,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.158262,'amu*angstrom^2'), symmetry=1, barrier=(6.31192,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.157649,'amu*angstrom^2'), symmetry=1, barrier=(6.31163,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.158478,'amu*angstrom^2'), symmetry=1, barrier=(6.31232,'kJ/mol'), semiclassical=False),
    ],
    spinMultiplicity = 2,
    opticalIsomers = 1,
    molecularWeight = (71.1408,'amu'),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[2.00959,0.0449637,-1.8456e-05,2.54053e-09,2.41293e-14,1470.95,19.172], Tmin=(100,'K'), Tmax=(1997.3,'K')), NASAPolynomial(coeffs=[14.1879,0.0277846,-1.09695e-05,1.84917e-09,-1.15577e-13,-4832.03,-51.6068], Tmin=(1997.3,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(11.6441,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(365.837,'J/(mol*K)'), comment="""Thermo group additivity estimation: group(Cs-CsCsCsH) + group(Cs-CsCsHH) + longDistanceInteraction_noncyclic(CsCs-ST) + group(Cs-CsHHH) + group(Cs-CsHHH) + group(Cs-CsHHH) + radical(Tertalkyl)"""),
)

species(
    label = 'C[C]C(114)',
    structure = SMILES('C[C]C'),
    E0 = (327.691,'kJ/mol'),
    modes = [
        HarmonicOscillator(frequencies=([2750,2770,2790,2810,2830,2850,1350,1400,1450,1500,700,800,1000,1100,1350,1400,900,1100,180],'cm^-1')),
        HinderedRotor(inertia=(0.0523925,'amu*angstrom^2'), symmetry=1, barrier=(10.4365,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.927751,'amu*angstrom^2'), symmetry=1, barrier=(21.3308,'kJ/mol'), semiclassical=False),
    ],
    spinMultiplicity = 3,
    opticalIsomers = 1,
    molecularWeight = (42.0797,'amu'),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[3.00517,0.0155111,1.85851e-05,-3.16862e-08,1.23732e-11,39453.6,10.1426], Tmin=(100,'K'), Tmax=(982.292,'K')), NASAPolynomial(coeffs=[6.73204,0.0159276,-5.86166e-06,1.06538e-09,-7.51285e-14,37969.2,-11.6], Tmin=(982.292,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(327.691,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(199.547,'J/(mol*K)'), comment="""Thermo group additivity estimation: group(Cs-CsCsHH) + group(Cs-CsHHH) + group(Cs-CsHHH) + radical(CCJ2_triplet)"""),
)

species(
    label = 'C[C]CCC(115)',
    structure = SMILES('C[C]CCC'),
    E0 = (280.131,'kJ/mol'),
    modes = [
        HarmonicOscillator(frequencies=([2750,2770,2790,2810,2830,2850,1350,1400,1450,1500,700,800,1000,1100,1350,1400,900,1100,2750,2783.33,2816.67,2850,1425,1450,1225,1275,1270,1340,700,800,300,400,274.483,274.534,3330.43],'cm^-1')),
        HinderedRotor(inertia=(0.117656,'amu*angstrom^2'), symmetry=1, barrier=(6.28725,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.454414,'amu*angstrom^2'), symmetry=1, barrier=(24.2856,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.454063,'amu*angstrom^2'), symmetry=1, barrier=(24.284,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.290904,'amu*angstrom^2'), symmetry=1, barrier=(15.5402,'kJ/mol'), semiclassical=False),
    ],
    spinMultiplicity = 3,
    opticalIsomers = 1,
    molecularWeight = (70.1329,'amu'),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[1.49126,0.0476291,-1.77616e-05,-4.75466e-09,3.6807e-12,33788.6,20.7635], Tmin=(100,'K'), Tmax=(1127.23,'K')), NASAPolynomial(coeffs=[10.1002,0.0302049,-1.20403e-05,2.19072e-09,-1.50455e-13,31013.9,-25.4993], Tmin=(1127.23,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(280.131,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(340.893,'J/(mol*K)'), comment="""Thermo group additivity estimation: group(Cs-CsCsHH) + group(Cs-CsCsHH) + group(Cs-CsCsHH) + group(Cs-CsHHH) + group(Cs-CsHHH) + radical(CCJ2_triplet)"""),
)

species(
    label = 'C=CCC(C)C(116)',
    structure = SMILES('C=CCC(C)C'),
    E0 = (-73.1238,'kJ/mol'),
    modes = [
        HarmonicOscillator(frequencies=([1380,1390,370,380,2900,435,2950,3100,1380,975,1025,1650,3010,987.5,1337.5,450,1655,2750,2850,1437.5,1250,1305,750,350,2750,2770,2790,2810,2830,2850,1350,1400,1450,1500,700,800,1000,1100,1350,1400,900,1100,180,4000],'cm^-1')),
        HinderedRotor(inertia=(0.830349,'amu*angstrom^2'), symmetry=1, barrier=(19.0914,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.830515,'amu*angstrom^2'), symmetry=1, barrier=(19.0952,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.0283688,'amu*angstrom^2'), symmetry=1, barrier=(12.3521,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.83097,'amu*angstrom^2'), symmetry=1, barrier=(19.1056,'kJ/mol'), semiclassical=False),
    ],
    spinMultiplicity = 1,
    opticalIsomers = 1,
    molecularWeight = (84.1595,'amu'),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[1.21741,0.0491153,7.0102e-06,-3.70579e-08,1.57922e-11,-8684.11,23.5832], Tmin=(100,'K'), Tmax=(1026.18,'K')), NASAPolynomial(coeffs=[11.4943,0.0363347,-1.41813e-05,2.61344e-09,-1.83421e-13,-12229.5,-33.2603], Tmin=(1026.18,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(-73.1238,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(415.724,'J/(mol*K)'), comment="""Thermo group additivity estimation: group(Cs-CsCsCsH) + group(Cs-(Cds-Cds)CsHH) + group(Cs-CsHHH) + group(Cs-CsHHH) + group(Cds-CdsCsH) + group(Cds-CdsHH)"""),
)

species(
    label = '[CH2]C(C)C(117)',
    structure = SMILES('[CH2]C(C)C'),
    E0 = (51.7368,'kJ/mol'),
    modes = [
        HarmonicOscillator(frequencies=([1380,1390,370,380,2900,435,3000,3100,440,815,1455,1000,2750,2770,2790,2810,2830,2850,1350,1400,1450,1500,700,800,1000,1100,1350,1400,900,1100],'cm^-1')),
        HinderedRotor(inertia=(0.670189,'amu*angstrom^2'), symmetry=1, barrier=(15.409,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.0486535,'amu*angstrom^2'), symmetry=1, barrier=(8.23719,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.011166,'amu*angstrom^2'), symmetry=1, barrier=(80.7214,'kJ/mol'), semiclassical=False),
    ],
    spinMultiplicity = 2,
    opticalIsomers = 1,
    molecularWeight = (57.1143,'amu'),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[2.40258,0.0268744,2.01407e-05,-4.094e-08,1.68979e-11,6287.47,16.146], Tmin=(100,'K'), Tmax=(953.315,'K')), NASAPolynomial(coeffs=[7.79053,0.024645,-8.4152e-06,1.45216e-09,-9.93173e-14,4334.2,-14.4467], Tmin=(953.315,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(51.7368,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(295.164,'J/(mol*K)'), comment="""Thermo group additivity estimation: group(Cs-CsCsCsH) + group(Cs-CsHHH) + group(Cs-CsHHH) + group(Cs-CsHHH) + radical(Isobutyl)"""),
)

species(
    label = 'C[CH]C(37)',
    structure = SMILES('C[CH]C'),
    E0 = (73.7998,'kJ/mol'),
    modes = [
        HarmonicOscillator(frequencies=([2750,2770,2790,2810,2830,2850,1350,1400,1450,1500,700,800,1000,1100,1350,1400,900,1100,3025,407.5,1350,352.5],'cm^-1')),
        HinderedRotor(inertia=(0.183496,'amu*angstrom^2'), symmetry=1, barrier=(4.21893,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.00321591,'amu*angstrom^2'), symmetry=1, barrier=(10.4862,'kJ/mol'), semiclassical=False),
    ],
    spinMultiplicity = 2,
    opticalIsomers = 1,
    molecularWeight = (43.0877,'amu'),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[3.21175,0.0131075,1.99174e-05,-2.55878e-08,8.25847e-12,8907.92,11.8466], Tmin=(100,'K'), Tmax=(1100.3,'K')), NASAPolynomial(coeffs=[3.65558,0.0230584,-9.41371e-06,1.73592e-09,-1.20102e-13,8110.22,6.48193], Tmin=(1100.3,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(73.7998,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(224.491,'J/(mol*K)'), comment="""Thermo group additivity estimation: group(Cs-CsCsHH) + group(Cs-CsHHH) + group(Cs-CsHHH) + radical(CCJC)"""),
)

species(
    label = '[CH2]C[CH]C(C)C(118)',
    structure = SMILES('[CH2]C[CH]C(C)C'),
    E0 = (202.23,'kJ/mol'),
    modes = [
        HarmonicOscillator(frequencies=([1380,1390,370,380,2900,435,3000,3100,440,815,1455,1000,2750,2850,1437.5,1250,1305,750,350,3025,407.5,1350,352.5,2750,2770,2790,2810,2830,2850,1350,1400,1450,1500,700,800,1000,1100,1350,1400,900,1100,1057.68,1058.04],'cm^-1')),
        HinderedRotor(inertia=(0.00379084,'amu*angstrom^2'), symmetry=1, barrier=(0.119627,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.167464,'amu*angstrom^2'), symmetry=1, barrier=(5.25426,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.00380077,'amu*angstrom^2'), symmetry=1, barrier=(0.119627,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.167503,'amu*angstrom^2'), symmetry=1, barrier=(5.26998,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.169019,'amu*angstrom^2'), symmetry=1, barrier=(5.2413,'kJ/mol'), semiclassical=False),
    ],
    spinMultiplicity = 3,
    opticalIsomers = 1,
    molecularWeight = (84.1595,'amu'),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[1.04857,0.0571479,-2.54632e-05,2.12507e-09,9.50513e-13,24435.2,28.0789], Tmin=(100,'K'), Tmax=(1373.69,'K')), NASAPolynomial(coeffs=[11.6436,0.0369441,-1.50284e-05,2.70351e-09,-1.81662e-13,20519.7,-30.0565], Tmin=(1373.69,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(202.23,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(411.566,'J/(mol*K)'), comment="""Thermo group additivity estimation: group(Cs-CsCsCsH) + group(Cs-CsCsHH) + longDistanceInteraction_noncyclic(CsCs-ST) + group(Cs-CsCsHH) + group(Cs-CsHHH) + group(Cs-CsHHH) + group(Cs-CsHHH) + radical(RCCJ) + radical(Cs_S)"""),
)

species(
    label = '[CH2][CH]CC(C)C(119)',
    structure = SMILES('[CH2][CH]CC(C)C'),
    E0 = (202.134,'kJ/mol'),
    modes = [
        HarmonicOscillator(frequencies=([1380,1390,370,380,2900,435,3000,3100,440,815,1455,1000,2750,2850,1437.5,1250,1305,750,350,3025,407.5,1350,352.5,2750,2770,2790,2810,2830,2850,1350,1400,1450,1500,700,800,1000,1100,1350,1400,900,1100,180,1079.34],'cm^-1')),
        HinderedRotor(inertia=(0.00608848,'amu*angstrom^2'), symmetry=1, barrier=(5.03337,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.218927,'amu*angstrom^2'), symmetry=1, barrier=(5.03356,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.0060888,'amu*angstrom^2'), symmetry=1, barrier=(5.03358,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.218926,'amu*angstrom^2'), symmetry=1, barrier=(5.03353,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.218917,'amu*angstrom^2'), symmetry=1, barrier=(5.03334,'kJ/mol'), semiclassical=False),
    ],
    spinMultiplicity = 3,
    opticalIsomers = 1,
    molecularWeight = (84.1595,'amu'),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[1.13689,0.0575366,-2.96173e-05,7.30425e-09,-7.19288e-13,24418.8,27.9492], Tmin=(100,'K'), Tmax=(2267.5,'K')), NASAPolynomial(coeffs=[17.2629,0.0290892,-1.07986e-05,1.77136e-09,-1.09264e-13,17105.7,-63.0514], Tmin=(2267.5,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(202.134,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(411.566,'J/(mol*K)'), comment="""Thermo group additivity estimation: group(Cs-CsCsCsH) + group(Cs-CsCsHH) + longDistanceInteraction_noncyclic(CsCs-ST) + group(Cs-CsCsHH) + group(Cs-CsHHH) + group(Cs-CsHHH) + group(Cs-CsHHH) + radical(RCCJ) + radical(RCCJC)"""),
)

species(
    label = '[CH2]CC(C)C(120)',
    structure = SMILES('[CH2]CC(C)C'),
    E0 = (31.4678,'kJ/mol'),
    modes = [
        HarmonicOscillator(frequencies=([1380,1390,370,380,2900,435,2750,2850,1437.5,1250,1305,750,350,3000,3100,440,815,1455,1000,2750,2770,2790,2810,2830,2850,1350,1400,1450,1500,700,800,1000,1100,1350,1400,900,1100,810.892],'cm^-1')),
        HinderedRotor(inertia=(0.132431,'amu*angstrom^2'), symmetry=1, barrier=(3.04485,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.47468,'amu*angstrom^2'), symmetry=1, barrier=(10.9138,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.132432,'amu*angstrom^2'), symmetry=1, barrier=(3.04488,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.0233793,'amu*angstrom^2'), symmetry=1, barrier=(10.9128,'kJ/mol'), semiclassical=False),
    ],
    spinMultiplicity = 2,
    opticalIsomers = 1,
    molecularWeight = (71.1408,'amu'),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[1.66349,0.0420054,3.89026e-06,-2.74708e-08,1.15723e-11,3876.87,21.0503], Tmin=(100,'K'), Tmax=(1042.13,'K')), NASAPolynomial(coeffs=[9.42897,0.0329343,-1.28984e-05,2.36148e-09,-1.64326e-13,1132.38,-22.1363], Tmin=(1042.13,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(31.4678,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(365.837,'J/(mol*K)'), comment="""Thermo group additivity estimation: group(Cs-CsCsCsH) + group(Cs-CsCsHH) + longDistanceInteraction_noncyclic(CsCs-ST) + group(Cs-CsHHH) + group(Cs-CsHHH) + group(Cs-CsHHH) + radical(RCCJ)"""),
)

species(
    label = '[CH]CCC(C)C(121)',
    structure = SMILES('[CH]CCC(C)C'),
    E0 = (250.656,'kJ/mol'),
    modes = [
        HarmonicOscillator(frequencies=([1380,1390,370,380,2900,435,2750,2783.33,2816.67,2850,1425,1450,1225,1275,1270,1340,700,800,300,400,2750,2770,2790,2810,2830,2850,1350,1400,1450,1500,700,800,1000,1100,1350,1400,900,1100,200,800,1066.67,1333.33,1600],'cm^-1')),
        HinderedRotor(inertia=(0.156089,'amu*angstrom^2'), symmetry=1, barrier=(3.5888,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.156089,'amu*angstrom^2'), symmetry=1, barrier=(3.5888,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.156089,'amu*angstrom^2'), symmetry=1, barrier=(3.5888,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.156089,'amu*angstrom^2'), symmetry=1, barrier=(3.5888,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.156089,'amu*angstrom^2'), symmetry=1, barrier=(3.5888,'kJ/mol'), semiclassical=False),
    ],
    spinMultiplicity = 3,
    opticalIsomers = 1,
    molecularWeight = (84.1595,'amu'),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[0.808226,0.0604888,-2.37075e-05,-6.10858e-09,5.08143e-12,30270.1,25.5846], Tmin=(100,'K'), Tmax=(1098.6,'K')), NASAPolynomial(coeffs=[12.2373,0.0362549,-1.43486e-05,2.6121e-09,-1.79974e-13,26710.1,-35.4024], Tmin=(1098.6,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(250.656,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(411.566,'J/(mol*K)'), comment="""Thermo group additivity estimation: group(Cs-CsCsCsH) + group(Cs-CsCsHH) + longDistanceInteraction_noncyclic(CsCs-ST) + group(Cs-CsCsHH) + group(Cs-CsHHH) + group(Cs-CsHHH) + group(Cs-CsHHH) + radical(CCJ2_triplet)"""),
)

species(
    label = 'CC=CC(C)C(122)',
    structure = SMILES('CC=CC(C)C'),
    E0 = (-84.2995,'kJ/mol'),
    modes = [
        HarmonicOscillator(frequencies=([1380,1390,370,380,2900,435,2995,3025,975,1000,1300,1375,400,500,1630,1680,2750,2762.5,2775,2787.5,2800,2812.5,2825,2837.5,2850,1350,1380,1410,1440,1470,1500,700,750,800,1000,1050,1100,1350,1375,1400,900,1000,1100,357.765],'cm^-1')),
        HinderedRotor(inertia=(0.113056,'amu*angstrom^2'), symmetry=1, barrier=(10.2776,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.113024,'amu*angstrom^2'), symmetry=1, barrier=(10.2778,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.113027,'amu*angstrom^2'), symmetry=1, barrier=(10.2793,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.113075,'amu*angstrom^2'), symmetry=1, barrier=(10.2801,'kJ/mol'), semiclassical=False),
    ],
    spinMultiplicity = 1,
    opticalIsomers = 1,
    molecularWeight = (84.1595,'amu'),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[1.26192,0.0490758,3.85704e-06,-3.18043e-08,1.3511e-11,-10030.7,23.0203], Tmin=(100,'K'), Tmax=(1042.84,'K')), NASAPolynomial(coeffs=[10.7305,0.0373741,-1.472e-05,2.70766e-09,-1.89056e-13,-13344.1,-29.4752], Tmin=(1042.84,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(-84.2995,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(415.724,'J/(mol*K)'), comment="""Thermo group additivity estimation: group(Cs-(Cds-Cds)CsCsH) + group(Cs-CsHHH) + group(Cs-CsHHH) + group(Cs-(Cds-Cds)HHH) + group(Cds-CdsCsH) + group(Cds-CdsCsH)"""),
)

species(
    label = 'CC=CCC(123)',
    structure = SMILES('CC=CCC'),
    E0 = (-52.5474,'kJ/mol'),
    modes = [
        HarmonicOscillator(frequencies=([2750,2770,2790,2810,2830,2850,1350,1400,1450,1500,700,800,1000,1100,1350,1400,900,1100,2750,2850,1437.5,1250,1305,750,350,2995,3025,975,1000,1300,1375,400,500,1630,1680,275.136],'cm^-1')),
        HinderedRotor(inertia=(0.170302,'amu*angstrom^2'), symmetry=1, barrier=(9.19048,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.170305,'amu*angstrom^2'), symmetry=1, barrier=(9.20606,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.172048,'amu*angstrom^2'), symmetry=1, barrier=(9.18149,'kJ/mol'), semiclassical=False),
    ],
    spinMultiplicity = 1,
    opticalIsomers = 1,
    molecularWeight = (70.1329,'amu'),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[1.94607,0.0364834,8.35486e-06,-2.83273e-08,1.11812e-11,-6238.66,19.2593], Tmin=(100,'K'), Tmax=(1059.79,'K')), NASAPolynomial(coeffs=[8.36587,0.0318147,-1.27245e-05,2.34965e-09,-1.63922e-13,-8697.93,-17.2682], Tmin=(1059.79,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(-52.5474,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(345.051,'J/(mol*K)'), comment="""Thermo group additivity estimation: group(Cs-(Cds-Cds)CsHH) + group(Cs-CsHHH) + group(Cs-(Cds-Cds)HHH) + group(Cds-CdsCsH) + group(Cds-CdsCsH)"""),
)

species(
    label = 'C=CC(C)C(124)',
    structure = SMILES('C=CC(C)C'),
    E0 = (-48.2739,'kJ/mol'),
    modes = [
        HarmonicOscillator(frequencies=([1380,1390,370,380,2900,435,2950,3100,1380,975,1025,1650,3010,987.5,1337.5,450,1655,2750,2770,2790,2810,2830,2850,1350,1400,1450,1500,700,800,1000,1100,1350,1400,900,1100,440.588],'cm^-1')),
        HinderedRotor(inertia=(0.0872468,'amu*angstrom^2'), symmetry=1, barrier=(12.0101,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.0867669,'amu*angstrom^2'), symmetry=1, barrier=(12.0158,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.0862256,'amu*angstrom^2'), symmetry=1, barrier=(12.0121,'kJ/mol'), semiclassical=False),
    ],
    spinMultiplicity = 1,
    opticalIsomers = 1,
    molecularWeight = (70.1329,'amu'),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[1.95972,0.0329298,2.62354e-05,-5.2288e-08,2.09885e-11,-5721.98,18.7076], Tmin=(100,'K'), Tmax=(989.921,'K')), NASAPolynomial(coeffs=[10.2073,0.0284297,-1.06251e-05,1.95192e-09,-1.38484e-13,-8767.28,-28.1327], Tmin=(989.921,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(-48.2739,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(345.051,'J/(mol*K)'), comment="""Thermo group additivity estimation: group(Cs-(Cds-Cds)CsCsH) + group(Cs-CsHHH) + group(Cs-CsHHH) + group(Cds-CdsCsH) + group(Cds-CdsHH)"""),
)

species(
    label = 'C[CH][CH]CC(125)',
    structure = SMILES('C[CH][CH]CC'),
    E0 = (220.82,'kJ/mol'),
    modes = [
        HarmonicOscillator(frequencies=([2750,2770,2790,2810,2830,2850,1350,1400,1450,1500,700,800,1000,1100,1350,1400,900,1100,2750,2850,1437.5,1250,1305,750,350,3000,3050,390,425,1340,1360,335,370,2185.7,2186.09],'cm^-1')),
        HinderedRotor(inertia=(0.0978821,'amu*angstrom^2'), symmetry=1, barrier=(2.2505,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.0966588,'amu*angstrom^2'), symmetry=1, barrier=(2.22238,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.0960164,'amu*angstrom^2'), symmetry=1, barrier=(2.2076,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.0966397,'amu*angstrom^2'), symmetry=1, barrier=(2.22194,'kJ/mol'), semiclassical=False),
    ],
    spinMultiplicity = 3,
    opticalIsomers = 1,
    molecularWeight = (70.1329,'amu'),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[2.82916,0.0356448,-8.78258e-06,-1.57857e-09,6.29749e-13,26591.6,21.3393], Tmin=(100,'K'), Tmax=(2058.92,'K')), NASAPolynomial(coeffs=[12.9595,0.0273799,-1.10783e-05,1.85779e-09,-1.14497e-13,20000.3,-40.7261], Tmin=(2058.92,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(220.82,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(340.893,'J/(mol*K)'), comment="""Thermo group additivity estimation: group(Cs-CsCsHH) + group(Cs-CsCsHH) + group(Cs-CsCsHH) + group(Cs-CsHHH) + group(Cs-CsHHH) + radical(RCCJCC) + radical(RCCJC)"""),
)

species(
    label = '[CH2][CH]C(C)C(126)',
    structure = SMILES('[CH2][CH]C(C)C'),
    E0 = (226.01,'kJ/mol'),
    modes = [
        HarmonicOscillator(frequencies=([1380,1390,370,380,2900,435,3025,407.5,1350,352.5,3000,3100,440,815,1455,1000,2750,2770,2790,2810,2830,2850,1350,1400,1450,1500,700,800,1000,1100,1350,1400,900,1100,1083.31],'cm^-1')),
        HinderedRotor(inertia=(0.290232,'amu*angstrom^2'), symmetry=1, barrier=(6.67301,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.144125,'amu*angstrom^2'), symmetry=1, barrier=(3.31371,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.136953,'amu*angstrom^2'), symmetry=1, barrier=(3.14882,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.00407018,'amu*angstrom^2'), symmetry=1, barrier=(3.37974,'kJ/mol'), semiclassical=False),
    ],
    spinMultiplicity = 3,
    opticalIsomers = 1,
    molecularWeight = (70.1329,'amu'),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[1.82817,0.0409722,-7.54625e-06,-1.01408e-08,4.45302e-12,27266.4,23.0232], Tmin=(100,'K'), Tmax=(1187.49,'K')), NASAPolynomial(coeffs=[8.34046,0.0322556,-1.32344e-05,2.42735e-09,-1.66588e-13,24787.7,-13.4385], Tmin=(1187.49,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(226.01,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(340.893,'J/(mol*K)'), comment="""Thermo group additivity estimation: group(Cs-CsCsCsH) + group(Cs-CsCsHH) + longDistanceInteraction_noncyclic(CsCs-ST) + group(Cs-CsHHH) + group(Cs-CsHHH) + group(Cs-CsHHH) + radical(RCCJ) + radical(Cs_S)"""),
)

species(
    label = 'C[CH][CH]C(C)C(127)',
    structure = SMILES('C[CH][CH]C(C)C'),
    E0 = (191.43,'kJ/mol'),
    modes = [
        HarmonicOscillator(frequencies=([1380,1390,370,380,2900,435,3000,3050,390,425,1340,1360,335,370,2750,2762.5,2775,2787.5,2800,2812.5,2825,2837.5,2850,1350,1380,1410,1440,1470,1500,700,750,800,1000,1050,1100,1350,1375,1400,900,1000,1100,209.994,1563.37],'cm^-1')),
        HinderedRotor(inertia=(0.118769,'amu*angstrom^2'), symmetry=1, barrier=(3.7166,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.118768,'amu*angstrom^2'), symmetry=1, barrier=(3.71661,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.118768,'amu*angstrom^2'), symmetry=1, barrier=(3.71661,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.118769,'amu*angstrom^2'), symmetry=1, barrier=(3.71661,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.00214286,'amu*angstrom^2'), symmetry=1, barrier=(3.71661,'kJ/mol'), semiclassical=False),
    ],
    spinMultiplicity = 3,
    opticalIsomers = 1,
    molecularWeight = (84.1595,'amu'),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[1.47072,0.0534612,-2.31224e-05,3.54482e-09,-1.18241e-14,23116.5,27.5206], Tmin=(100,'K'), Tmax=(1857.96,'K')), NASAPolynomial(coeffs=[14.3789,0.0334697,-1.32786e-05,2.2718e-09,-1.44505e-13,16973.9,-46.3725], Tmin=(1857.96,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(191.43,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(411.566,'J/(mol*K)'), comment="""Thermo group additivity estimation: group(Cs-CsCsCsH) + group(Cs-CsCsHH) + longDistanceInteraction_noncyclic(CsCs-ST) + group(Cs-CsCsHH) + group(Cs-CsHHH) + group(Cs-CsHHH) + group(Cs-CsHHH) + radical(Cs_S) + radical(RCCJC)"""),
)

species(
    label = 'C[CH]C(C)C(128)',
    structure = SMILES('C[CH]C(C)C'),
    E0 = (20.7636,'kJ/mol'),
    modes = [
        HarmonicOscillator(frequencies=([1380,1390,370,380,2900,435,3025,407.5,1350,352.5,2750,2762.5,2775,2787.5,2800,2812.5,2825,2837.5,2850,1350,1380,1410,1440,1470,1500,700,750,800,1000,1050,1100,1350,1375,1400,900,1000,1100,180.005],'cm^-1')),
        HinderedRotor(inertia=(0.0716601,'amu*angstrom^2'), symmetry=1, barrier=(1.64761,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.0716832,'amu*angstrom^2'), symmetry=1, barrier=(1.64857,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.0683943,'amu*angstrom^2'), symmetry=1, barrier=(1.57252,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.0699751,'amu*angstrom^2'), symmetry=1, barrier=(1.60895,'kJ/mol'), semiclassical=False),
    ],
    spinMultiplicity = 2,
    opticalIsomers = 1,
    molecularWeight = (71.1408,'amu'),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[1.80209,0.040159,2.93649e-06,-2.20734e-08,8.5979e-12,2583.19,21.3269], Tmin=(100,'K'), Tmax=(1112.67,'K')), NASAPolynomial(coeffs=[8.12229,0.0352151,-1.4364e-05,2.6515e-09,-1.8373e-13,76.3184,-14.7841], Tmin=(1112.67,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(20.7636,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(365.837,'J/(mol*K)'), comment="""Thermo group additivity estimation: group(Cs-CsCsCsH) + group(Cs-CsCsHH) + longDistanceInteraction_noncyclic(CsCs-ST) + group(Cs-CsHHH) + group(Cs-CsHHH) + group(Cs-CsHHH) + radical(Cs_S)"""),
)

species(
    label = '[CH]C(C)C(131)',
    structure = SMILES('[CH]C(C)C'),
    E0 = (294.87,'kJ/mol'),
    modes = [
        HarmonicOscillator(frequencies=([1380,1390,370,380,2900,435,2750,2770,2790,2810,2830,2850,1350,1400,1450,1500,700,800,1000,1100,1350,1400,900,1100,180,1069.26,2660.96],'cm^-1')),
        HinderedRotor(inertia=(0.770701,'amu*angstrom^2'), symmetry=1, barrier=(17.7199,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.00285457,'amu*angstrom^2'), symmetry=1, barrier=(2.35755,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.770149,'amu*angstrom^2'), symmetry=1, barrier=(17.7072,'kJ/mol'), semiclassical=False),
    ],
    spinMultiplicity = 3,
    opticalIsomers = 1,
    molecularWeight = (56.1063,'amu'),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[2.32395,0.0283397,1.27929e-05,-3.33016e-08,1.39108e-11,35532.5,15.6509], Tmin=(100,'K'), Tmax=(989.36,'K')), NASAPolynomial(coeffs=[8.92727,0.0218887,-8.12276e-06,1.47626e-09,-1.03817e-13,33235,-21.1434], Tmin=(989.36,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(294.87,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(270.22,'J/(mol*K)'), comment="""Thermo group additivity estimation: group(Cs-CsCsCsH) + group(Cs-CsHHH) + group(Cs-CsHHH) + group(Cs-CsHHH) + radical(CCJ2_triplet)"""),
)

species(
    label = 'CC[C]C(C)C(132)',
    structure = SMILES('CC[C]C(C)C'),
    E0 = (250.656,'kJ/mol'),
    modes = [
        HarmonicOscillator(frequencies=([1380,1390,370,380,2900,435,2750,2850,1437.5,1250,1305,750,350,2750,2762.5,2775,2787.5,2800,2812.5,2825,2837.5,2850,1350,1380,1410,1440,1470,1500,700,750,800,1000,1050,1100,1350,1375,1400,900,1000,1100,200,800,1600],'cm^-1')),
        HinderedRotor(inertia=(0.156089,'amu*angstrom^2'), symmetry=1, barrier=(3.5888,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.156089,'amu*angstrom^2'), symmetry=1, barrier=(3.5888,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.156089,'amu*angstrom^2'), symmetry=1, barrier=(3.5888,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.156089,'amu*angstrom^2'), symmetry=1, barrier=(3.5888,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.156089,'amu*angstrom^2'), symmetry=1, barrier=(3.5888,'kJ/mol'), semiclassical=False),
    ],
    spinMultiplicity = 3,
    opticalIsomers = 1,
    molecularWeight = (84.1595,'amu'),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[0.808226,0.0604888,-2.37075e-05,-6.10858e-09,5.08143e-12,30270.1,24.4859], Tmin=(100,'K'), Tmax=(1098.6,'K')), NASAPolynomial(coeffs=[12.2373,0.0362549,-1.43486e-05,2.6121e-09,-1.79974e-13,26710.1,-36.501], Tmin=(1098.6,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(250.656,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(411.566,'J/(mol*K)'), comment="""Thermo group additivity estimation: group(Cs-CsCsCsH) + group(Cs-CsCsHH) + longDistanceInteraction_noncyclic(CsCs-ST) + group(Cs-CsCsHH) + group(Cs-CsHHH) + group(Cs-CsHHH) + group(Cs-CsHHH) + radical(CCJ2_triplet)"""),
)

species(
    label = 'C=C(CC)CC(133)',
    structure = SMILES('C=C(CC)CC'),
    E0 = (-78.2223,'kJ/mol'),
    modes = [
        HarmonicOscillator(frequencies=([2750,2770,2790,2810,2830,2850,1350,1400,1450,1500,700,800,1000,1100,1350,1400,900,1100,2950,3100,1380,975,1025,1650,2750,2783.33,2816.67,2850,1425,1450,1225,1275,1270,1340,700,800,300,400,350,440,435,1725,357.457,357.458],'cm^-1')),
        HinderedRotor(inertia=(0.00131935,'amu*angstrom^2'), symmetry=1, barrier=(0.119627,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.119871,'amu*angstrom^2'), symmetry=1, barrier=(10.8687,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.119867,'amu*angstrom^2'), symmetry=1, barrier=(10.8687,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.119867,'amu*angstrom^2'), symmetry=1, barrier=(10.8687,'kJ/mol'), semiclassical=False),
    ],
    spinMultiplicity = 1,
    opticalIsomers = 1,
    molecularWeight = (84.1595,'amu'),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[1.10231,0.05276,-3.99671e-06,-2.53431e-08,1.15942e-11,-9294.24,23.4572], Tmin=(100,'K'), Tmax=(1047.64,'K')), NASAPolynomial(coeffs=[11.4876,0.0365644,-1.43928e-05,2.6442e-09,-1.84401e-13,-12757.5,-33.2727], Tmin=(1047.64,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(-78.2223,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(415.724,'J/(mol*K)'), comment="""Thermo group additivity estimation: group(Cs-(Cds-Cds)CsHH) + group(Cs-(Cds-Cds)CsHH) + group(Cs-CsHHH) + group(Cs-CsHHH) + group(Cds-CdsCsCs) + group(Cds-CdsHH)"""),
)

species(
    label = '[CH2][C](CC)CC(136)',
    structure = SMILES('[CH2][C](CC)CC'),
    E0 = (196.293,'kJ/mol'),
    modes = [
        HarmonicOscillator(frequencies=([2750,2770,2790,2810,2830,2850,1350,1400,1450,1500,700,800,1000,1100,1350,1400,900,1100,2750,2783.33,2816.67,2850,1425,1450,1225,1275,1270,1340,700,800,300,400,3000,3100,440,815,1455,1000,360,370,350,180,2476.06],'cm^-1')),
        HinderedRotor(inertia=(0.144223,'amu*angstrom^2'), symmetry=1, barrier=(3.31598,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.144225,'amu*angstrom^2'), symmetry=1, barrier=(3.31601,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.144238,'amu*angstrom^2'), symmetry=1, barrier=(3.31631,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.144204,'amu*angstrom^2'), symmetry=1, barrier=(3.31554,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.144245,'amu*angstrom^2'), symmetry=1, barrier=(3.31648,'kJ/mol'), semiclassical=False),
    ],
    spinMultiplicity = 3,
    opticalIsomers = 1,
    molecularWeight = (84.1595,'amu'),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[2.03401,0.0509563,1.07262e-05,-9.44335e-08,8.02874e-11,23671.6,23.2587], Tmin=(100,'K'), Tmax=(471.198,'K')), NASAPolynomial(coeffs=[3.66017,0.0490308,-2.09591e-05,3.89767e-09,-2.69542e-13,23386.4,15.2377], Tmin=(471.198,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(196.293,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(411.566,'J/(mol*K)'), comment="""Thermo group additivity estimation: group(Cs-CsCsCsH) + group(Cs-CsCsHH) + longDistanceInteraction_noncyclic(CsCs-ST) + group(Cs-CsCsHH) + longDistanceInteraction_noncyclic(CsCs-ST) + group(Cs-CsHHH) + group(Cs-CsHHH) + group(Cs-CsHHH) + radical(Isobutyl) + radical(Tertalkyl)"""),
)

species(
    label = '[CH2]C([CH2])CC(137)',
    structure = SMILES('[CH2]C([CH2])CC'),
    E0 = (236.386,'kJ/mol'),
    modes = [
        HarmonicOscillator(frequencies=([1380,1390,370,380,2900,435,2750,2850,1437.5,1250,1305,750,350,3000,3033.33,3066.67,3100,415,465,780,850,1435,1475,900,1100,2750,2800,2850,1350,1500,750,1050,1375,1000,1687.38],'cm^-1')),
        HinderedRotor(inertia=(0.076419,'amu*angstrom^2'), symmetry=1, barrier=(8.22007,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.00407207,'amu*angstrom^2'), symmetry=1, barrier=(8.22277,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.0766863,'amu*angstrom^2'), symmetry=1, barrier=(8.22257,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.647884,'amu*angstrom^2'), symmetry=1, barrier=(69.6769,'kJ/mol'), semiclassical=False),
    ],
    spinMultiplicity = 3,
    opticalIsomers = 1,
    molecularWeight = (70.1329,'amu'),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[1.71506,0.0440446,-8.68986e-06,-1.71203e-08,9.74755e-12,28518.6,22.6405], Tmin=(100,'K'), Tmax=(931.108,'K')), NASAPolynomial(coeffs=[8.79894,0.0291761,-9.80972e-06,1.63338e-09,-1.07806e-13,26524.8,-14.6523], Tmin=(931.108,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(236.386,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(340.893,'J/(mol*K)'), comment="""Thermo group additivity estimation: group(Cs-CsCsCsH) + group(Cs-CsCsHH) + longDistanceInteraction_noncyclic(CsCs-ST) + group(Cs-CsHHH) + group(Cs-CsHHH) + group(Cs-CsHHH) + radical(Isobutyl) + radical(Isobutyl)"""),
)

species(
    label = '[CH2]C([CH]C)CC(138)',
    structure = SMILES('[CH2]C([CH]C)CC'),
    E0 = (205.413,'kJ/mol'),
    modes = [
        HarmonicOscillator(frequencies=([1380,1390,370,380,2900,435,3000,3100,440,815,1455,1000,2750,2850,1437.5,1250,1305,750,350,3025,407.5,1350,352.5,2750,2770,2790,2810,2830,2850,1350,1400,1450,1500,700,800,1000,1100,1350,1400,900,1100,180,1249.07],'cm^-1')),
        HinderedRotor(inertia=(0.0904743,'amu*angstrom^2'), symmetry=1, barrier=(2.08018,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.0912646,'amu*angstrom^2'), symmetry=1, barrier=(2.09835,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.0900749,'amu*angstrom^2'), symmetry=1, barrier=(2.071,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.0908681,'amu*angstrom^2'), symmetry=1, barrier=(2.08924,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.0904328,'amu*angstrom^2'), symmetry=1, barrier=(2.07923,'kJ/mol'), semiclassical=False),
    ],
    spinMultiplicity = 3,
    opticalIsomers = 1,
    molecularWeight = (84.1595,'amu'),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[0.867061,0.0598704,-3.31987e-05,9.27409e-09,-1.05724e-12,24826,29.4302], Tmin=(100,'K'), Tmax=(1969.66,'K')), NASAPolynomial(coeffs=[15.0165,0.0311353,-1.13151e-05,1.8671e-09,-1.17089e-13,19252.2,-48.4237], Tmin=(1969.66,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(205.413,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(411.566,'J/(mol*K)'), comment="""Thermo group additivity estimation: group(Cs-CsCsCsH) + group(Cs-CsCsHH) + longDistanceInteraction_noncyclic(CsCs-ST) + group(Cs-CsCsHH) + longDistanceInteraction_noncyclic(CsCs-ST) + group(Cs-CsHHH) + group(Cs-CsHHH) + group(Cs-CsHHH) + radical(Isobutyl) + radical(Cs_S)"""),
)

species(
    label = '[CH2]CC([CH2])CC(139)',
    structure = SMILES('[CH2]CC([CH2])CC'),
    E0 = (216.117,'kJ/mol'),
    modes = [
        HarmonicOscillator(frequencies=([1380,1390,370,380,2900,435,2750,2783.33,2816.67,2850,1425,1450,1225,1275,1270,1340,700,800,300,400,3000,3033.33,3066.67,3100,415,465,780,850,1435,1475,900,1100,2750,2800,2850,1350,1500,750,1050,1375,1000,1075.22,2799.79],'cm^-1')),
        HinderedRotor(inertia=(0.10897,'amu*angstrom^2'), symmetry=1, barrier=(2.50544,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.107115,'amu*angstrom^2'), symmetry=1, barrier=(2.46279,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.691222,'amu*angstrom^2'), symmetry=1, barrier=(15.8926,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.108477,'amu*angstrom^2'), symmetry=1, barrier=(2.4941,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(3.51947,'amu*angstrom^2'), symmetry=1, barrier=(80.9194,'kJ/mol'), semiclassical=False),
    ],
    spinMultiplicity = 3,
    opticalIsomers = 1,
    molecularWeight = (84.1595,'amu'),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[0.905828,0.0600018,-2.7815e-05,3.76906e-11,2.86971e-12,26111,28.4894], Tmin=(100,'K'), Tmax=(1094.78,'K')), NASAPolynomial(coeffs=[10.8091,0.0368474,-1.39422e-05,2.46082e-09,-1.66083e-13,23161.8,-23.7512], Tmin=(1094.78,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(216.117,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(411.566,'J/(mol*K)'), comment="""Thermo group additivity estimation: group(Cs-CsCsCsH) + group(Cs-CsCsHH) + longDistanceInteraction_noncyclic(CsCs-ST) + group(Cs-CsCsHH) + longDistanceInteraction_noncyclic(CsCs-ST) + group(Cs-CsHHH) + group(Cs-CsHHH) + group(Cs-CsHHH) + radical(Isobutyl) + radical(RCCJ)"""),
)

species(
    label = '[CH]C(CC)CC(140)',
    structure = SMILES('[CH]C(CC)CC'),
    E0 = (254.004,'kJ/mol'),
    modes = [
        HarmonicOscillator(frequencies=([1380,1390,370,380,2900,435,2750,2783.33,2816.67,2850,1425,1450,1225,1275,1270,1340,700,800,300,400,2750,2770,2790,2810,2830,2850,1350,1400,1450,1500,700,800,1000,1100,1350,1400,900,1100,200,800,1066.67,1333.33,1600],'cm^-1')),
        HinderedRotor(inertia=(0.156089,'amu*angstrom^2'), symmetry=1, barrier=(3.5888,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.156089,'amu*angstrom^2'), symmetry=1, barrier=(3.5888,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.156089,'amu*angstrom^2'), symmetry=1, barrier=(3.5888,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.156089,'amu*angstrom^2'), symmetry=1, barrier=(3.5888,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.156089,'amu*angstrom^2'), symmetry=1, barrier=(3.5888,'kJ/mol'), semiclassical=False),
    ],
    spinMultiplicity = 3,
    opticalIsomers = 1,
    molecularWeight = (84.1595,'amu'),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[0.808226,0.0604888,-2.37075e-05,-6.10858e-09,5.08143e-12,30672.6,25.5846], Tmin=(100,'K'), Tmax=(1098.6,'K')), NASAPolynomial(coeffs=[12.2373,0.0362549,-1.43486e-05,2.6121e-09,-1.79974e-13,27112.7,-35.4024], Tmin=(1098.6,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(254.004,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(411.566,'J/(mol*K)'), comment="""Thermo group additivity estimation: group(Cs-CsCsCsH) + group(Cs-CsCsHH) + longDistanceInteraction_noncyclic(CsCs-ST) + group(Cs-CsCsHH) + longDistanceInteraction_noncyclic(CsCs-ST) + group(Cs-CsHHH) + group(Cs-CsHHH) + group(Cs-CsHHH) + radical(CCJ2_triplet)"""),
)

species(
    label = 'C[CH]C[CH]C(141)',
    structure = SMILES('C[CH]C[CH]C'),
    E0 = (220.808,'kJ/mol'),
    modes = [
        HarmonicOscillator(frequencies=([2750,2770,2790,2810,2830,2850,1350,1400,1450,1500,700,800,1000,1100,1350,1400,900,1100,2750,2850,1437.5,1250,1305,750,350,3000,3050,390,425,1340,1360,335,370,1975.88,1978.06],'cm^-1')),
        HinderedRotor(inertia=(0.119273,'amu*angstrom^2'), symmetry=1, barrier=(2.74233,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.118618,'amu*angstrom^2'), symmetry=1, barrier=(2.72726,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.119225,'amu*angstrom^2'), symmetry=1, barrier=(2.74122,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.117483,'amu*angstrom^2'), symmetry=1, barrier=(2.70116,'kJ/mol'), semiclassical=False),
    ],
    spinMultiplicity = 3,
    opticalIsomers = 1,
    molecularWeight = (70.1329,'amu'),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[2.5835,0.0375164,-1.14014e-05,-2.10758e-10,4.03354e-13,26602.2,21.7084], Tmin=(100,'K'), Tmax=(2050.38,'K')), NASAPolynomial(coeffs=[12.6275,0.02724,-1.07003e-05,1.7777e-09,-1.09341e-13,20524.7,-38.7365], Tmin=(2050.38,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(220.808,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(340.893,'J/(mol*K)'), comment="""Thermo group additivity estimation: group(Cs-CsCsHH) + group(Cs-CsCsHH) + group(Cs-CsCsHH) + group(Cs-CsHHH) + group(Cs-CsHHH) + radical(RCCJC) + radical(RCCJC)"""),
)

species(
    label = '[CH]CC(C)C(142)',
    structure = SMILES('[CH]CC(C)C'),
    E0 = (274.437,'kJ/mol'),
    modes = [
        HarmonicOscillator(frequencies=([1380,1390,370,380,2900,435,2750,2850,1437.5,1250,1305,750,350,2750,2770,2790,2810,2830,2850,1350,1400,1450,1500,700,800,1000,1100,1350,1400,900,1100,180,1178.08,1181.01,2850.55],'cm^-1')),
        HinderedRotor(inertia=(0.20908,'amu*angstrom^2'), symmetry=1, barrier=(4.80716,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.721932,'amu*angstrom^2'), symmetry=1, barrier=(16.5986,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.0167809,'amu*angstrom^2'), symmetry=1, barrier=(16.5991,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.720321,'amu*angstrom^2'), symmetry=1, barrier=(16.5616,'kJ/mol'), semiclassical=False),
    ],
    spinMultiplicity = 3,
    opticalIsomers = 1,
    molecularWeight = (70.1329,'amu'),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[1.57675,0.0443051,-5.15955e-06,-1.99594e-08,9.54331e-12,33102.1,20.5784], Tmin=(100,'K'), Tmax=(1033.91,'K')), NASAPolynomial(coeffs=[10.4344,0.0293104,-1.13681e-05,2.07458e-09,-1.44363e-13,30240.3,-27.432], Tmin=(1033.91,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(274.437,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(340.893,'J/(mol*K)'), comment="""Thermo group additivity estimation: group(Cs-CsCsCsH) + group(Cs-CsCsHH) + longDistanceInteraction_noncyclic(CsCs-ST) + group(Cs-CsHHH) + group(Cs-CsHHH) + group(Cs-CsHHH) + radical(CCJ2_triplet)"""),
)

species(
    label = 'C[C]CC(C)C(143)',
    structure = SMILES('C[C]CC(C)C'),
    E0 = (250.656,'kJ/mol'),
    modes = [
        HarmonicOscillator(frequencies=([1380,1390,370,380,2900,435,2750,2850,1437.5,1250,1305,750,350,2750,2762.5,2775,2787.5,2800,2812.5,2825,2837.5,2850,1350,1380,1410,1440,1470,1500,700,750,800,1000,1050,1100,1350,1375,1400,900,1000,1100,200,800,1600],'cm^-1')),
        HinderedRotor(inertia=(0.156089,'amu*angstrom^2'), symmetry=1, barrier=(3.5888,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.156089,'amu*angstrom^2'), symmetry=1, barrier=(3.5888,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.156089,'amu*angstrom^2'), symmetry=1, barrier=(3.5888,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.156089,'amu*angstrom^2'), symmetry=1, barrier=(3.5888,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.156089,'amu*angstrom^2'), symmetry=1, barrier=(3.5888,'kJ/mol'), semiclassical=False),
    ],
    spinMultiplicity = 3,
    opticalIsomers = 1,
    molecularWeight = (84.1595,'amu'),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[0.808226,0.0604888,-2.37075e-05,-6.10858e-09,5.08143e-12,30270.1,24.4859], Tmin=(100,'K'), Tmax=(1098.6,'K')), NASAPolynomial(coeffs=[12.2373,0.0362549,-1.43486e-05,2.6121e-09,-1.79974e-13,26710.1,-36.501], Tmin=(1098.6,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(250.656,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(411.566,'J/(mol*K)'), comment="""Thermo group additivity estimation: group(Cs-CsCsCsH) + group(Cs-CsCsHH) + longDistanceInteraction_noncyclic(CsCs-ST) + group(Cs-CsCsHH) + group(Cs-CsHHH) + group(Cs-CsHHH) + group(Cs-CsHHH) + radical(CCJ2_triplet)"""),
)

species(
    label = 'CC=C(C)CC(144)',
    structure = SMILES('CC=C(C)CC'),
    E0 = (-91.6024,'kJ/mol'),
    modes = [
        HarmonicOscillator(frequencies=([2750,2762.5,2775,2787.5,2800,2812.5,2825,2837.5,2850,1350,1380,1410,1440,1470,1500,700,750,800,1000,1050,1100,1350,1375,1400,900,1000,1100,3010,987.5,1337.5,450,1655,2750,2850,1437.5,1250,1305,750,350,350,440,435,1725,262.84],'cm^-1')),
        HinderedRotor(inertia=(0.174652,'amu*angstrom^2'), symmetry=1, barrier=(8.48771,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.177628,'amu*angstrom^2'), symmetry=1, barrier=(8.51283,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.174547,'amu*angstrom^2'), symmetry=1, barrier=(8.50099,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.177543,'amu*angstrom^2'), symmetry=1, barrier=(8.50725,'kJ/mol'), semiclassical=False),
    ],
    spinMultiplicity = 1,
    opticalIsomers = 1,
    molecularWeight = (84.1595,'amu'),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[1.02076,0.0572653,-2.40349e-05,-1.96769e-10,1.90923e-12,-10903.3,23.3489], Tmin=(100,'K'), Tmax=(1260.33,'K')), NASAPolynomial(coeffs=[11.0038,0.0376957,-1.51618e-05,2.73618e-09,-1.85344e-13,-14381.8,-30.9404], Tmin=(1260.33,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(-91.6024,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(415.724,'J/(mol*K)'), comment="""Thermo group additivity estimation: group(Cs-(Cds-Cds)CsHH) + group(Cs-CsHHH) + group(Cs-(Cds-Cds)HHH) + group(Cs-(Cds-Cds)HHH) + group(Cds-CdsCsCs) + group(Cds-CdsCsH)"""),
)

species(
    label = 'C=C(C)CC(145)',
    structure = SMILES('C=C(C)CC'),
    E0 = (-55.5768,'kJ/mol'),
    modes = [
        HarmonicOscillator(frequencies=([2750,2770,2790,2810,2830,2850,1350,1400,1450,1500,700,800,1000,1100,1350,1400,900,1100,2950,3100,1380,975,1025,1650,2750,2850,1437.5,1250,1305,750,350,350,440,435,1725,368.476],'cm^-1')),
        HinderedRotor(inertia=(0.100804,'amu*angstrom^2'), symmetry=1, barrier=(9.71422,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.100839,'amu*angstrom^2'), symmetry=1, barrier=(9.71429,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.10097,'amu*angstrom^2'), symmetry=1, barrier=(9.71361,'kJ/mol'), semiclassical=False),
    ],
    spinMultiplicity = 1,
    opticalIsomers = 1,
    molecularWeight = (70.1329,'amu'),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[1.77093,0.0406537,-7.47548e-07,-2.08242e-08,9.00712e-12,-6596.98,18.8383], Tmin=(100,'K'), Tmax=(1062.08,'K')), NASAPolynomial(coeffs=[9.1127,0.0308933,-1.22293e-05,2.24281e-09,-1.55787e-13,-9165.5,-21.7739], Tmin=(1062.08,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(-55.5768,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(345.051,'J/(mol*K)'), comment="""Thermo group additivity estimation: group(Cs-(Cds-Cds)CsHH) + group(Cs-CsHHH) + group(Cs-(Cds-Cds)HHH) + group(Cds-CdsCsCs) + group(Cds-CdsHH)"""),
)

species(
    label = '[CH2][C](C)CC(146)',
    structure = SMILES('[CH2][C](C)CC'),
    E0 = (216.726,'kJ/mol'),
    modes = [
        HarmonicOscillator(frequencies=([2750,2770,2790,2810,2830,2850,1350,1400,1450,1500,700,800,1000,1100,1350,1400,900,1100,2750,2850,1437.5,1250,1305,750,350,3000,3100,440,815,1455,1000,360,370,350,1818.46],'cm^-1')),
        HinderedRotor(inertia=(0.0027399,'amu*angstrom^2'), symmetry=1, barrier=(6.42812,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.279805,'amu*angstrom^2'), symmetry=1, barrier=(6.43326,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.279899,'amu*angstrom^2'), symmetry=1, barrier=(6.43543,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.2798,'amu*angstrom^2'), symmetry=1, barrier=(6.43316,'kJ/mol'), semiclassical=False),
    ],
    spinMultiplicity = 3,
    opticalIsomers = 1,
    molecularWeight = (70.1329,'amu'),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[1.94841,0.0477752,-3.58852e-05,2.3051e-08,-7.78947e-12,26137.5,21.8541], Tmin=(100,'K'), Tmax=(777.642,'K')), NASAPolynomial(coeffs=[3.0281,0.0400455,-1.6778e-05,3.07209e-09,-2.09742e-13,26035.4,17.3398], Tmin=(777.642,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(216.726,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(340.893,'J/(mol*K)'), comment="""Thermo group additivity estimation: group(Cs-CsCsCsH) + group(Cs-CsCsHH) + longDistanceInteraction_noncyclic(CsCs-ST) + group(Cs-CsHHH) + group(Cs-CsHHH) + group(Cs-CsHHH) + radical(Tertalkyl) + radical(Isobutyl)"""),
)

species(
    label = 'C[CH][C](C)CC(147)',
    structure = SMILES('C[CH][C](C)CC'),
    E0 = (185.753,'kJ/mol'),
    modes = [
        HarmonicOscillator(frequencies=([2750,2762.5,2775,2787.5,2800,2812.5,2825,2837.5,2850,1350,1380,1410,1440,1470,1500,700,750,800,1000,1050,1100,1350,1375,1400,900,1000,1100,2750,2850,1437.5,1250,1305,750,350,3025,407.5,1350,352.5,360,370,350,180,3774.43],'cm^-1')),
        HinderedRotor(inertia=(0.0961369,'amu*angstrom^2'), symmetry=1, barrier=(2.21038,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.09613,'amu*angstrom^2'), symmetry=1, barrier=(2.21022,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.096127,'amu*angstrom^2'), symmetry=1, barrier=(2.21015,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.0961443,'amu*angstrom^2'), symmetry=1, barrier=(2.21055,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.0961466,'amu*angstrom^2'), symmetry=1, barrier=(2.2106,'kJ/mol'), semiclassical=False),
    ],
    spinMultiplicity = 3,
    opticalIsomers = 1,
    molecularWeight = (84.1595,'amu'),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[2.19765,0.0516106,-2.20518e-05,3.48251e-09,-1.22645e-13,22394.1,23.9393], Tmin=(100,'K'), Tmax=(2362.29,'K')), NASAPolynomial(coeffs=[25.0431,0.0221863,-9.24759e-06,1.52827e-09,-9.1428e-14,9017.05,-111.384], Tmin=(2362.29,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(185.753,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(411.566,'J/(mol*K)'), comment="""Thermo group additivity estimation: group(Cs-CsCsCsH) + group(Cs-CsCsHH) + longDistanceInteraction_noncyclic(CsCs-ST) + group(Cs-CsCsHH) + longDistanceInteraction_noncyclic(CsCs-ST) + group(Cs-CsHHH) + group(Cs-CsHHH) + group(Cs-CsHHH) + radical(Tertalkyl) + radical(Cs_S)"""),
)

species(
    label = '[CH2]C[C](C)CC(148)',
    structure = SMILES('[CH2]C[C](C)CC'),
    E0 = (196.457,'kJ/mol'),
    modes = [
        HarmonicOscillator(frequencies=([2750,2770,2790,2810,2830,2850,1350,1400,1450,1500,700,800,1000,1100,1350,1400,900,1100,2750,2783.33,2816.67,2850,1425,1450,1225,1275,1270,1340,700,800,300,400,3000,3100,440,815,1455,1000,360,370,350,180,2549.89],'cm^-1')),
        HinderedRotor(inertia=(0.133245,'amu*angstrom^2'), symmetry=1, barrier=(3.06356,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.133245,'amu*angstrom^2'), symmetry=1, barrier=(3.06357,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.13324,'amu*angstrom^2'), symmetry=1, barrier=(3.06344,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.13324,'amu*angstrom^2'), symmetry=1, barrier=(3.06345,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.133244,'amu*angstrom^2'), symmetry=1, barrier=(3.06353,'kJ/mol'), semiclassical=False),
    ],
    spinMultiplicity = 3,
    opticalIsomers = 1,
    molecularWeight = (84.1595,'amu'),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[1.97979,0.0547295,-2.65263e-05,5.76402e-09,-4.86714e-13,23689.2,23.9148], Tmin=(100,'K'), Tmax=(2565.56,'K')), NASAPolynomial(coeffs=[17.1226,0.03112,-1.27225e-05,2.17706e-09,-1.37183e-13,15919.3,-63.4079], Tmin=(2565.56,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(196.457,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(411.566,'J/(mol*K)'), comment="""Thermo group additivity estimation: group(Cs-CsCsCsH) + group(Cs-CsCsHH) + longDistanceInteraction_noncyclic(CsCs-ST) + group(Cs-CsCsHH) + longDistanceInteraction_noncyclic(CsCs-ST) + group(Cs-CsHHH) + group(Cs-CsHHH) + group(Cs-CsHHH) + radical(Tertalkyl) + radical(RCCJ)"""),
)

species(
    label = 'C[C]CC(149)',
    structure = SMILES('C[C]CC'),
    E0 = (303.911,'kJ/mol'),
    modes = [
        HarmonicOscillator(frequencies=([2750,2770,2790,2810,2830,2850,1350,1400,1450,1500,700,800,1000,1100,1350,1400,900,1100,2750,2850,1437.5,1250,1305,750,350,180,976.109],'cm^-1')),
        HinderedRotor(inertia=(0.727183,'amu*angstrom^2'), symmetry=1, barrier=(16.7194,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.726898,'amu*angstrom^2'), symmetry=1, barrier=(16.7128,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.220844,'amu*angstrom^2'), symmetry=1, barrier=(5.07763,'kJ/mol'), semiclassical=False),
    ],
    spinMultiplicity = 3,
    opticalIsomers = 1,
    molecularWeight = (56.1063,'amu'),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[2.2609,0.0314424,7.48735e-07,-1.84873e-08,8.06281e-12,36620.5,15.7527], Tmin=(100,'K'), Tmax=(1038.41,'K')), NASAPolynomial(coeffs=[8.22231,0.0233776,-9.12316e-06,1.66749e-09,-1.15986e-13,34579.2,-17.1001], Tmin=(1038.41,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(303.911,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(270.22,'J/(mol*K)'), comment="""Thermo group additivity estimation: group(Cs-CsCsHH) + group(Cs-CsCsHH) + group(Cs-CsHHH) + group(Cs-CsHHH) + radical(CCJ2_triplet)"""),
)

species(
    label = 'CC[C]CC(150)',
    structure = SMILES('CC[C]CC'),
    E0 = (280.131,'kJ/mol'),
    modes = [
        HarmonicOscillator(frequencies=([2750,2770,2790,2810,2830,2850,1350,1400,1450,1500,700,800,1000,1100,1350,1400,900,1100,2750,2783.33,2816.67,2850,1425,1450,1225,1275,1270,1340,700,800,300,400,274.483,274.534,3330.43],'cm^-1')),
        HinderedRotor(inertia=(0.117656,'amu*angstrom^2'), symmetry=1, barrier=(6.28725,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.454414,'amu*angstrom^2'), symmetry=1, barrier=(24.2856,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.454063,'amu*angstrom^2'), symmetry=1, barrier=(24.284,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.290904,'amu*angstrom^2'), symmetry=1, barrier=(15.5402,'kJ/mol'), semiclassical=False),
    ],
    spinMultiplicity = 3,
    opticalIsomers = 1,
    molecularWeight = (70.1329,'amu'),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[1.49126,0.0476291,-1.77616e-05,-4.75466e-09,3.6807e-12,33788.6,20.0704], Tmin=(100,'K'), Tmax=(1127.23,'K')), NASAPolynomial(coeffs=[10.1002,0.0302049,-1.20403e-05,2.19072e-09,-1.50455e-13,31013.9,-26.1924], Tmin=(1127.23,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(280.131,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(340.893,'J/(mol*K)'), comment="""Thermo group additivity estimation: group(Cs-CsCsHH) + group(Cs-CsCsHH) + group(Cs-CsCsHH) + group(Cs-CsHHH) + group(Cs-CsHHH) + radical(CCJ2_triplet)"""),
)

species(
    label = 'C=CC(C)CC(151)',
    structure = SMILES('C=CC(C)CC'),
    E0 = (-72.0542,'kJ/mol'),
    modes = [
        HarmonicOscillator(frequencies=([1380,1390,370,380,2900,435,2950,3100,1380,975,1025,1650,3010,987.5,1337.5,450,1655,2750,2850,1437.5,1250,1305,750,350,2750,2770,2790,2810,2830,2850,1350,1400,1450,1500,700,800,1000,1100,1350,1400,900,1100,225.831,225.837],'cm^-1')),
        HinderedRotor(inertia=(0.00330569,'amu*angstrom^2'), symmetry=1, barrier=(0.119627,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.00330587,'amu*angstrom^2'), symmetry=1, barrier=(0.119627,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.401592,'amu*angstrom^2'), symmetry=1, barrier=(14.5349,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.401685,'amu*angstrom^2'), symmetry=1, barrier=(14.5349,'kJ/mol'), semiclassical=False),
    ],
    spinMultiplicity = 1,
    opticalIsomers = 1,
    molecularWeight = (84.1595,'amu'),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[1.21824,0.0488307,8.49214e-06,-3.91884e-08,1.67101e-11,-8555.22,24.3076], Tmin=(100,'K'), Tmax=(1018.92,'K')), NASAPolynomial(coeffs=[11.6699,0.0359249,-1.39119e-05,2.55989e-09,-1.7982e-13,-12145.1,-33.4759], Tmin=(1018.92,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(-72.0542,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(415.724,'J/(mol*K)'), comment="""Thermo group additivity estimation: group(Cs-(Cds-Cds)CsCsH) + group(Cs-CsCsHH) + group(Cs-CsHHH) + group(Cs-CsHHH) + group(Cds-CdsCsH) + group(Cds-CdsHH)"""),
)

species(
    label = 'CC=CC(152)',
    structure = SMILES('CC=CC'),
    E0 = (-29.9019,'kJ/mol'),
    modes = [
        HarmonicOscillator(frequencies=([2750,2770,2790,2810,2830,2850,1350,1400,1450,1500,700,800,1000,1100,1350,1400,900,1100,2995,3025,975,1000,1300,1375,400,500,1630,1680],'cm^-1')),
        HinderedRotor(inertia=(0.470151,'amu*angstrom^2'), symmetry=1, barrier=(10.8097,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.472334,'amu*angstrom^2'), symmetry=1, barrier=(10.8599,'kJ/mol'), semiclassical=False),
    ],
    spinMultiplicity = 1,
    opticalIsomers = 1,
    molecularWeight = (56.1063,'amu'),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[2.61745,0.0243422,1.17366e-05,-2.39951e-08,8.67979e-12,-3541.51,13.2443], Tmin=(100,'K'), Tmax=(1081.87,'K')), NASAPolynomial(coeffs=[5.99499,0.0261386,-1.05588e-05,1.94786e-09,-1.35282e-13,-5108.26,-7.17958], Tmin=(1081.87,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(-29.9019,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(274.378,'J/(mol*K)'), comment="""Thermo group additivity estimation: group(Cs-(Cds-Cds)HHH) + group(Cs-(Cds-Cds)HHH) + group(Cds-CdsCsH) + group(Cds-CdsCsH)"""),
)

species(
    label = 'C[CH][CH]C(153)',
    structure = SMILES('C[CH][CH]C'),
    E0 = (244.589,'kJ/mol'),
    modes = [
        HarmonicOscillator(frequencies=([2750,2770,2790,2810,2830,2850,1350,1400,1450,1500,700,800,1000,1100,1350,1400,900,1100,3000,3050,390,425,1340,1360,335,370,180],'cm^-1')),
        HinderedRotor(inertia=(0.00322957,'amu*angstrom^2'), symmetry=1, barrier=(6.49009,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.00322929,'amu*angstrom^2'), symmetry=1, barrier=(6.49394,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.00322879,'amu*angstrom^2'), symmetry=1, barrier=(6.49067,'kJ/mol'), semiclassical=False),
    ],
    spinMultiplicity = 3,
    opticalIsomers = 1,
    molecularWeight = (56.1063,'amu'),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[3.1236,0.0238473,-9.09934e-07,-4.51998e-09,1.14171e-12,29444.7,17.5352], Tmin=(100,'K'), Tmax=(1813.48,'K')), NASAPolynomial(coeffs=[6.29329,0.026079,-1.03847e-05,1.76763e-09,-1.11706e-13,26778.4,-3.825], Tmin=(1813.48,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(244.589,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(270.22,'J/(mol*K)'), comment="""Thermo group additivity estimation: group(Cs-CsCsHH) + group(Cs-CsCsHH) + group(Cs-CsHHH) + group(Cs-CsHHH) + radical(RCCJC) + radical(RCCJC)"""),
)

species(
    label = '[CH2]C(C)[CH]C(154)',
    structure = SMILES('[CH2]C(C)[CH]C'),
    E0 = (225.846,'kJ/mol'),
    modes = [
        HarmonicOscillator(frequencies=([1380,1390,370,380,2900,435,3025,407.5,1350,352.5,3000,3100,440,815,1455,1000,2750,2770,2790,2810,2830,2850,1350,1400,1450,1500,700,800,1000,1100,1350,1400,900,1100,1290.55],'cm^-1')),
        HinderedRotor(inertia=(0.0998922,'amu*angstrom^2'), symmetry=1, barrier=(2.29672,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.0999483,'amu*angstrom^2'), symmetry=1, barrier=(2.29801,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.099759,'amu*angstrom^2'), symmetry=1, barrier=(2.29366,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.0998509,'amu*angstrom^2'), symmetry=1, barrier=(2.29577,'kJ/mol'), semiclassical=False),
    ],
    spinMultiplicity = 3,
    opticalIsomers = 1,
    molecularWeight = (70.1329,'amu'),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[1.80046,0.0421184,-1.07296e-05,-7.76563e-09,4.08541e-12,27247.4,23.8047], Tmin=(100,'K'), Tmax=(1119,'K')), NASAPolynomial(coeffs=[7.75177,0.0320962,-1.23774e-05,2.20164e-09,-1.48914e-13,25211,-8.72364], Tmin=(1119,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(225.846,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(340.893,'J/(mol*K)'), comment="""Thermo group additivity estimation: group(Cs-CsCsCsH) + group(Cs-CsCsHH) + longDistanceInteraction_noncyclic(CsCs-ST) + group(Cs-CsHHH) + group(Cs-CsHHH) + group(Cs-CsHHH) + radical(Isobutyl) + radical(Cs_S)"""),
)

species(
    label = 'C[CH]C(C)[CH]C(155)',
    structure = SMILES('C[CH]C(C)[CH]C'),
    E0 = (194.873,'kJ/mol'),
    modes = [
        HarmonicOscillator(frequencies=([1380,1390,370,380,2900,435,3000,3050,390,425,1340,1360,335,370,2750,2762.5,2775,2787.5,2800,2812.5,2825,2837.5,2850,1350,1380,1410,1440,1470,1500,700,750,800,1000,1050,1100,1350,1375,1400,900,1000,1100,1159.51,4000],'cm^-1')),
        HinderedRotor(inertia=(0.207382,'amu*angstrom^2'), symmetry=1, barrier=(5.98228,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.207385,'amu*angstrom^2'), symmetry=1, barrier=(5.98227,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.00627035,'amu*angstrom^2'), symmetry=1, barrier=(5.98227,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(1.04309,'amu*angstrom^2'), symmetry=1, barrier=(30.0894,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.00414704,'amu*angstrom^2'), symmetry=1, barrier=(0.119627,'kJ/mol'), semiclassical=False),
    ],
    spinMultiplicity = 3,
    opticalIsomers = 1,
    molecularWeight = (84.1595,'amu'),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[1.29964,0.0539078,-2.13068e-05,6.97715e-10,9.11612e-13,23539.4,27.9576], Tmin=(100,'K'), Tmax=(1530.47,'K')), NASAPolynomial(coeffs=[11.9179,0.0370216,-1.54057e-05,2.76585e-09,-1.83922e-13,19016.7,-31.9453], Tmin=(1530.47,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(194.873,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(411.566,'J/(mol*K)'), comment="""Thermo group additivity estimation: group(Cs-CsCsCsH) + group(Cs-CsCsHH) + longDistanceInteraction_noncyclic(CsCs-ST) + group(Cs-CsCsHH) + longDistanceInteraction_noncyclic(CsCs-ST) + group(Cs-CsHHH) + group(Cs-CsHHH) + group(Cs-CsHHH) + radical(Cs_S) + radical(Cs_S)"""),
)

species(
    label = '[CH2]CC(C)[CH]C(156)',
    structure = SMILES('[CH2]CC(C)[CH]C'),
    E0 = (205.577,'kJ/mol'),
    modes = [
        HarmonicOscillator(frequencies=([1380,1390,370,380,2900,435,3000,3100,440,815,1455,1000,2750,2850,1437.5,1250,1305,750,350,3025,407.5,1350,352.5,2750,2770,2790,2810,2830,2850,1350,1400,1450,1500,700,800,1000,1100,1350,1400,900,1100,1057.79,1057.91],'cm^-1')),
        HinderedRotor(inertia=(0.00381087,'amu*angstrom^2'), symmetry=1, barrier=(0.119627,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.167176,'amu*angstrom^2'), symmetry=1, barrier=(5.2554,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.00381922,'amu*angstrom^2'), symmetry=1, barrier=(0.119627,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.167385,'amu*angstrom^2'), symmetry=1, barrier=(5.25483,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.167776,'amu*angstrom^2'), symmetry=1, barrier=(5.25489,'kJ/mol'), semiclassical=False),
    ],
    spinMultiplicity = 3,
    opticalIsomers = 1,
    molecularWeight = (84.1595,'amu'),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[1.04855,0.057148,-2.54636e-05,2.12545e-09,9.50389e-13,24837.7,28.7721], Tmin=(100,'K'), Tmax=(1373.7,'K')), NASAPolynomial(coeffs=[11.6438,0.0369438,-1.50283e-05,2.70347e-09,-1.81659e-13,20922.2,-29.3645], Tmin=(1373.7,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(205.577,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(411.566,'J/(mol*K)'), comment="""Thermo group additivity estimation: group(Cs-CsCsCsH) + group(Cs-CsCsHH) + longDistanceInteraction_noncyclic(CsCs-ST) + group(Cs-CsCsHH) + longDistanceInteraction_noncyclic(CsCs-ST) + group(Cs-CsHHH) + group(Cs-CsHHH) + group(Cs-CsHHH) + radical(Cs_S) + radical(RCCJ)"""),
)

species(
    label = '[CH2][CH]C(C)CC(157)',
    structure = SMILES('[CH2][CH]C(C)CC'),
    E0 = (205.577,'kJ/mol'),
    modes = [
        HarmonicOscillator(frequencies=([1380,1390,370,380,2900,435,3000,3100,440,815,1455,1000,2750,2850,1437.5,1250,1305,750,350,3025,407.5,1350,352.5,2750,2770,2790,2810,2830,2850,1350,1400,1450,1500,700,800,1000,1100,1350,1400,900,1100,1057.79,1057.91],'cm^-1')),
        HinderedRotor(inertia=(0.00381087,'amu*angstrom^2'), symmetry=1, barrier=(0.119627,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.167176,'amu*angstrom^2'), symmetry=1, barrier=(5.2554,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.00381922,'amu*angstrom^2'), symmetry=1, barrier=(0.119627,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.167385,'amu*angstrom^2'), symmetry=1, barrier=(5.25483,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.167776,'amu*angstrom^2'), symmetry=1, barrier=(5.25489,'kJ/mol'), semiclassical=False),
    ],
    spinMultiplicity = 3,
    opticalIsomers = 1,
    molecularWeight = (84.1595,'amu'),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[1.04855,0.057148,-2.54636e-05,2.12545e-09,9.50389e-13,24837.7,28.7721], Tmin=(100,'K'), Tmax=(1373.7,'K')), NASAPolynomial(coeffs=[11.6438,0.0369438,-1.50283e-05,2.70347e-09,-1.81659e-13,20922.2,-29.3645], Tmin=(1373.7,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(205.577,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(411.566,'J/(mol*K)'), comment="""Thermo group additivity estimation: group(Cs-CsCsCsH) + group(Cs-CsCsHH) + longDistanceInteraction_noncyclic(CsCs-ST) + group(Cs-CsCsHH) + longDistanceInteraction_noncyclic(CsCs-ST) + group(Cs-CsHHH) + group(Cs-CsHHH) + group(Cs-CsHHH) + radical(Cs_S) + radical(RCCJ)"""),
)

species(
    label = 'C[CH]CC(43)',
    structure = SMILES('C[CH]CC'),
    E0 = (50.1423,'kJ/mol'),
    modes = [
        HarmonicOscillator(frequencies=([2750,2770,2790,2810,2830,2850,1350,1400,1450,1500,700,800,1000,1100,1350,1400,900,1100,2750,2850,1437.5,1250,1305,750,350,3025,407.5,1350,352.5,1960.71],'cm^-1')),
        HinderedRotor(inertia=(0.180902,'amu*angstrom^2'), symmetry=1, barrier=(4.1593,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.181078,'amu*angstrom^2'), symmetry=1, barrier=(4.16334,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.180901,'amu*angstrom^2'), symmetry=1, barrier=(4.15927,'kJ/mol'), semiclassical=False),
    ],
    spinMultiplicity = 2,
    opticalIsomers = 1,
    molecularWeight = (57.1143,'amu'),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[2.51758,0.0284121,2.01392e-06,-1.19574e-08,4.04516e-12,6087.4,17.6691], Tmin=(100,'K'), Tmax=(1230.11,'K')), NASAPolynomial(coeffs=[5.03064,0.0300836,-1.20272e-05,2.15735e-09,-1.45453e-13,4724.39,1.99745], Tmin=(1230.11,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(50.1423,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(295.164,'J/(mol*K)'), comment="""Thermo group additivity estimation: group(Cs-CsCsHH) + group(Cs-CsCsHH) + group(Cs-CsHHH) + group(Cs-CsHHH) + radical(RCCJC)"""),
)

species(
    label = '[CH]C(C)CC(158)',
    structure = SMILES('[CH]C(C)CC'),
    E0 = (274.437,'kJ/mol'),
    modes = [
        HarmonicOscillator(frequencies=([1380,1390,370,380,2900,435,2750,2850,1437.5,1250,1305,750,350,2750,2770,2790,2810,2830,2850,1350,1400,1450,1500,700,800,1000,1100,1350,1400,900,1100,180,1178.08,1181.01,2850.55],'cm^-1')),
        HinderedRotor(inertia=(0.20908,'amu*angstrom^2'), symmetry=1, barrier=(4.80716,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.721932,'amu*angstrom^2'), symmetry=1, barrier=(16.5986,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.0167809,'amu*angstrom^2'), symmetry=1, barrier=(16.5991,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.720321,'amu*angstrom^2'), symmetry=1, barrier=(16.5616,'kJ/mol'), semiclassical=False),
    ],
    spinMultiplicity = 3,
    opticalIsomers = 1,
    molecularWeight = (70.1329,'amu'),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[1.57675,0.0443051,-5.15955e-06,-1.99594e-08,9.54331e-12,33102.1,21.2716], Tmin=(100,'K'), Tmax=(1033.91,'K')), NASAPolynomial(coeffs=[10.4344,0.0293104,-1.13681e-05,2.07458e-09,-1.44363e-13,30240.3,-26.7388], Tmin=(1033.91,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(274.437,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(340.893,'J/(mol*K)'), comment="""Thermo group additivity estimation: group(Cs-CsCsCsH) + group(Cs-CsCsHH) + longDistanceInteraction_noncyclic(CsCs-ST) + group(Cs-CsHHH) + group(Cs-CsHHH) + group(Cs-CsHHH) + radical(CCJ2_triplet)"""),
)

species(
    label = 'C[C]C(C)CC(159)',
    structure = SMILES('C[C]C(C)CC'),
    E0 = (254.004,'kJ/mol'),
    modes = [
        HarmonicOscillator(frequencies=([1380,1390,370,380,2900,435,2750,2850,1437.5,1250,1305,750,350,2750,2762.5,2775,2787.5,2800,2812.5,2825,2837.5,2850,1350,1380,1410,1440,1470,1500,700,750,800,1000,1050,1100,1350,1375,1400,900,1000,1100,200,800,1600],'cm^-1')),
        HinderedRotor(inertia=(0.156089,'amu*angstrom^2'), symmetry=1, barrier=(3.5888,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.156089,'amu*angstrom^2'), symmetry=1, barrier=(3.5888,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.156089,'amu*angstrom^2'), symmetry=1, barrier=(3.5888,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.156089,'amu*angstrom^2'), symmetry=1, barrier=(3.5888,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.156089,'amu*angstrom^2'), symmetry=1, barrier=(3.5888,'kJ/mol'), semiclassical=False),
    ],
    spinMultiplicity = 3,
    opticalIsomers = 1,
    molecularWeight = (84.1595,'amu'),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[0.808226,0.0604888,-2.37075e-05,-6.10858e-09,5.08143e-12,30672.6,25.1791], Tmin=(100,'K'), Tmax=(1098.6,'K')), NASAPolynomial(coeffs=[12.2373,0.0362549,-1.43486e-05,2.6121e-09,-1.79974e-13,27112.7,-35.8078], Tmin=(1098.6,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(254.004,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(411.566,'J/(mol*K)'), comment="""Thermo group additivity estimation: group(Cs-CsCsCsH) + group(Cs-CsCsHH) + longDistanceInteraction_noncyclic(CsCs-ST) + group(Cs-CsCsHH) + longDistanceInteraction_noncyclic(CsCs-ST) + group(Cs-CsHHH) + group(Cs-CsHHH) + group(Cs-CsHHH) + radical(CCJ2_triplet)"""),
)

species(
    label = '[CH2]C([CH2])(C)C(162)',
    structure = SMILES('[CH2]C([CH2])(C)C'),
    E0 = (220.81,'kJ/mol'),
    modes = [
        HarmonicOscillator(frequencies=([350,400,1190,1240,450,3000,3033.33,3066.67,3100,415,465,780,850,1435,1475,900,1100,2750,2770,2790,2810,2830,2850,1350,1400,1450,1500,700,800,1000,1100,1350,1400,900,1100],'cm^-1')),
        HinderedRotor(inertia=(0.566762,'amu*angstrom^2'), symmetry=1, barrier=(13.031,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.0755856,'amu*angstrom^2'), symmetry=1, barrier=(1.73786,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.0735248,'amu*angstrom^2'), symmetry=1, barrier=(1.69048,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.566599,'amu*angstrom^2'), symmetry=1, barrier=(13.0272,'kJ/mol'), semiclassical=False),
    ],
    spinMultiplicity = 3,
    opticalIsomers = 1,
    molecularWeight = (70.1329,'amu'),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[1.50184,0.0464035,-1.23397e-05,-1.22214e-08,6.69061e-12,26654.6,20.5795], Tmin=(100,'K'), Tmax=(1068.51,'K')), NASAPolynomial(coeffs=[10.7678,0.0287282,-1.14087e-05,2.09835e-09,-1.46078e-13,23703.3,-29.2819], Tmin=(1068.51,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(220.81,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(340.893,'J/(mol*K)'), comment="""Thermo group additivity estimation: group(Cs-CsCsCsCs) + group(Cs-CsHHH) + group(Cs-CsHHH) + group(Cs-CsHHH) + group(Cs-CsHHH) + radical(Neopentyl) + radical(Neopentyl)"""),
)

species(
    label = '[CH2]C(C)(C)[CH]C(163)',
    structure = SMILES('[CH2]C(C)(C)[CH]C'),
    E0 = (193.269,'kJ/mol'),
    modes = [
        HarmonicOscillator(frequencies=([350,400,1190,1240,450,3025,407.5,1350,352.5,3000,3100,440,815,1455,1000,2750,2762.5,2775,2787.5,2800,2812.5,2825,2837.5,2850,1350,1380,1410,1440,1470,1500,700,750,800,1000,1050,1100,1350,1375,1400,900,1000,1100,760.042],'cm^-1')),
        HinderedRotor(inertia=(0.114038,'amu*angstrom^2'), symmetry=1, barrier=(2.62196,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.118719,'amu*angstrom^2'), symmetry=1, barrier=(2.72959,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.112553,'amu*angstrom^2'), symmetry=1, barrier=(2.58782,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.54045,'amu*angstrom^2'), symmetry=1, barrier=(12.426,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.537538,'amu*angstrom^2'), symmetry=1, barrier=(12.3591,'kJ/mol'), semiclassical=False),
    ],
    spinMultiplicity = 3,
    opticalIsomers = 1,
    molecularWeight = (84.1595,'amu'),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[0.921091,0.0583337,-2.24476e-05,-4.13839e-09,3.52761e-12,23363.4,26.4542], Tmin=(100,'K'), Tmax=(1201.55,'K')), NASAPolynomial(coeffs=[12.1012,0.0368787,-1.53429e-05,2.83855e-09,-1.95875e-13,19538.7,-34.2715], Tmin=(1201.55,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(193.269,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(411.566,'J/(mol*K)'), comment="""Thermo group additivity estimation: group(Cs-CsCsCsCs) + group(Cs-CsCsHH) + longDistanceInteraction_noncyclic(CsCs-SQ) + group(Cs-CsHHH) + group(Cs-CsHHH) + group(Cs-CsHHH) + group(Cs-CsHHH) + radical(Cs_S) + radical(Neopentyl)"""),
)

species(
    label = '[CH2]C([CH2])(C)CC(164)',
    structure = SMILES('[CH2]C([CH2])(C)CC'),
    E0 = (203.725,'kJ/mol'),
    modes = [
        HarmonicOscillator(frequencies=([350,400,1190,1240,450,2750,2850,1437.5,1250,1305,750,350,3000,3033.33,3066.67,3100,415,465,780,850,1435,1475,900,1100,2750,2770,2790,2810,2830,2850,1350,1400,1450,1500,700,800,1000,1100,1350,1400,900,1100,4000],'cm^-1')),
        HinderedRotor(inertia=(0.18176,'amu*angstrom^2'), symmetry=1, barrier=(4.17903,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.182087,'amu*angstrom^2'), symmetry=1, barrier=(4.18654,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.984864,'amu*angstrom^2'), symmetry=1, barrier=(22.644,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.671713,'amu*angstrom^2'), symmetry=1, barrier=(15.444,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.985811,'amu*angstrom^2'), symmetry=1, barrier=(22.6657,'kJ/mol'), semiclassical=False),
    ],
    spinMultiplicity = 3,
    opticalIsomers = 1,
    molecularWeight = (84.1595,'amu'),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[0.717573,0.0627359,-3.12312e-05,1.82149e-09,2.25434e-12,24628.4,26.3377], Tmin=(100,'K'), Tmax=(1163.44,'K')), NASAPolynomial(coeffs=[12.8968,0.0351592,-1.41095e-05,2.57254e-09,-1.76606e-13,20826.8,-38.4222], Tmin=(1163.44,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(203.725,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(411.566,'J/(mol*K)'), comment="""Thermo group additivity estimation: group(Cs-CsCsCsCs) + group(Cs-CsCsHH) + longDistanceInteraction_noncyclic(CsCs-SQ) + group(Cs-CsHHH) + group(Cs-CsHHH) + group(Cs-CsHHH) + group(Cs-CsHHH) + radical(Neopentyl) + radical(Neopentyl)"""),
)

species(
    label = '[CH2]CC([CH2])(C)C(165)',
    structure = SMILES('[CH2]CC([CH2])(C)C'),
    E0 = (203.973,'kJ/mol'),
    modes = [
        HarmonicOscillator(frequencies=([350,400,1190,1240,450,2750,2850,1437.5,1250,1305,750,350,3000,3033.33,3066.67,3100,415,465,780,850,1435,1475,900,1100,2750,2770,2790,2810,2830,2850,1350,1400,1450,1500,700,800,1000,1100,1350,1400,900,1100,180],'cm^-1')),
        HinderedRotor(inertia=(0.115583,'amu*angstrom^2'), symmetry=1, barrier=(2.65749,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.0141314,'amu*angstrom^2'), symmetry=1, barrier=(2.64806,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.115737,'amu*angstrom^2'), symmetry=1, barrier=(2.66101,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.692476,'amu*angstrom^2'), symmetry=1, barrier=(15.9214,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.692389,'amu*angstrom^2'), symmetry=1, barrier=(15.9194,'kJ/mol'), semiclassical=False),
    ],
    spinMultiplicity = 3,
    opticalIsomers = 1,
    molecularWeight = (84.1595,'amu'),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[0.790922,0.0601355,-2.15902e-05,-9.06187e-09,6.15828e-12,24656.6,26.1436], Tmin=(100,'K'), Tmax=(1097.58,'K')), NASAPolynomial(coeffs=[12.9329,0.0353335,-1.42729e-05,2.6372e-09,-1.83534e-13,20819.8,-38.9011], Tmin=(1097.58,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(203.973,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(411.566,'J/(mol*K)'), comment="""Thermo group additivity estimation: group(Cs-CsCsCsCs) + group(Cs-CsCsHH) + longDistanceInteraction_noncyclic(CsCs-SQ) + group(Cs-CsHHH) + group(Cs-CsHHH) + group(Cs-CsHHH) + group(Cs-CsHHH) + radical(Neopentyl) + radical(RCCJ)"""),
)

species(
    label = '[CH2]C(C)(C)C(166)',
    structure = SMILES('[CH2]C(C)(C)C'),
    E0 = (15.8124,'kJ/mol'),
    modes = [
        HarmonicOscillator(frequencies=([350,400,1190,1240,450,3000,3100,440,815,1455,1000,2750,2762.5,2775,2787.5,2800,2812.5,2825,2837.5,2850,1350,1380,1410,1440,1470,1500,700,750,800,1000,1050,1100,1350,1375,1400,900,1000,1100],'cm^-1')),
        HinderedRotor(inertia=(0.535044,'amu*angstrom^2'), symmetry=1, barrier=(12.3017,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.534071,'amu*angstrom^2'), symmetry=1, barrier=(12.2793,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.28524,'amu*angstrom^2'), symmetry=1, barrier=(12.2391,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.0330443,'amu*angstrom^2'), symmetry=1, barrier=(12.1928,'kJ/mol'), semiclassical=False),
    ],
    spinMultiplicity = 2,
    opticalIsomers = 1,
    molecularWeight = (71.1408,'amu'),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[1.51889,0.0432662,7.19675e-06,-3.48195e-08,1.48765e-11,2000.93,18.3974], Tmin=(100,'K'), Tmax=(1023.92,'K')), NASAPolynomial(coeffs=[11.3066,0.0307292,-1.20856e-05,2.24784e-09,-1.58974e-13,-1350.6,-35.6325], Tmin=(1023.92,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(15.8124,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(365.837,'J/(mol*K)'), comment="""Thermo group additivity estimation: group(Cs-CsCsCsCs) + group(Cs-CsHHH) + group(Cs-CsHHH) + group(Cs-CsHHH) + group(Cs-CsHHH) + radical(Neopentyl)"""),
)

species(
    label = '[CH]C(C)(C)CC(167)',
    structure = SMILES('[CH]C(C)(C)CC'),
    E0 = (241.944,'kJ/mol'),
    modes = [
        HarmonicOscillator(frequencies=([350,400,1190,1240,450,2750,2850,1437.5,1250,1305,750,350,2750,2762.5,2775,2787.5,2800,2812.5,2825,2837.5,2850,1350,1380,1410,1440,1470,1500,700,750,800,1000,1050,1100,1350,1375,1400,900,1000,1100,200,800,1200,1600],'cm^-1')),
        HinderedRotor(inertia=(0.156089,'amu*angstrom^2'), symmetry=1, barrier=(3.5888,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.156089,'amu*angstrom^2'), symmetry=1, barrier=(3.5888,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.156089,'amu*angstrom^2'), symmetry=1, barrier=(3.5888,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.156089,'amu*angstrom^2'), symmetry=1, barrier=(3.5888,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.156089,'amu*angstrom^2'), symmetry=1, barrier=(3.5888,'kJ/mol'), semiclassical=False),
    ],
    spinMultiplicity = 3,
    opticalIsomers = 1,
    molecularWeight = (84.1595,'amu'),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[0.728645,0.0592265,-1.09334e-05,-2.42486e-08,1.23035e-11,29227.8,23.8675], Tmin=(100,'K'), Tmax=(1039.07,'K')), NASAPolynomial(coeffs=[14.341,0.033928,-1.35393e-05,2.52714e-09,-1.78676e-13,24935.9,-49.3667], Tmin=(1039.07,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(241.944,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(411.566,'J/(mol*K)'), comment="""Thermo group additivity estimation: group(Cs-CsCsCsCs) + group(Cs-CsCsHH) + longDistanceInteraction_noncyclic(CsCs-SQ) + group(Cs-CsHHH) + group(Cs-CsHHH) + group(Cs-CsHHH) + group(Cs-CsHHH) + radical(CCJ2_triplet)"""),
)

species(
    label = '[CH2]CC(C)C[CH2](168)',
    structure = SMILES('[CH2]CC(C)C[CH2]'),
    E0 = (216.281,'kJ/mol'),
    modes = [
        HarmonicOscillator(frequencies=([1380,1390,370,380,2900,435,2750,2783.33,2816.67,2850,1425,1450,1225,1275,1270,1340,700,800,300,400,3000,3033.33,3066.67,3100,415,465,780,850,1435,1475,900,1100,2750,2800,2850,1350,1500,750,1050,1375,1000,180,1163.55],'cm^-1')),
        HinderedRotor(inertia=(0.140192,'amu*angstrom^2'), symmetry=1, barrier=(3.22328,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.138487,'amu*angstrom^2'), symmetry=1, barrier=(3.1841,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.0033347,'amu*angstrom^2'), symmetry=1, barrier=(3.19479,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.138453,'amu*angstrom^2'), symmetry=1, barrier=(3.1833,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.777352,'amu*angstrom^2'), symmetry=1, barrier=(17.8729,'kJ/mol'), semiclassical=False),
    ],
    spinMultiplicity = 3,
    opticalIsomers = 1,
    molecularWeight = (84.1595,'amu'),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[0.908203,0.059148,-2.56192e-05,-1.11615e-09,2.74554e-12,26131.2,27.7991], Tmin=(100,'K'), Tmax=(1174.87,'K')), NASAPolynomial(coeffs=[11.5146,0.0368117,-1.46882e-05,2.66059e-09,-1.81623e-13,22688.3,-29.1258], Tmin=(1174.87,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(216.281,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(411.566,'J/(mol*K)'), comment="""Thermo group additivity estimation: group(Cs-CsCsCsH) + group(Cs-CsCsHH) + longDistanceInteraction_noncyclic(CsCs-ST) + group(Cs-CsCsHH) + longDistanceInteraction_noncyclic(CsCs-ST) + group(Cs-CsHHH) + group(Cs-CsHHH) + group(Cs-CsHHH) + radical(RCCJ) + radical(RCCJ)"""),
)

species(
    label = '[CH]CC(C)CC(169)',
    structure = SMILES('[CH]CC(C)CC'),
    E0 = (254.004,'kJ/mol'),
    modes = [
        HarmonicOscillator(frequencies=([1380,1390,370,380,2900,435,2750,2783.33,2816.67,2850,1425,1450,1225,1275,1270,1340,700,800,300,400,2750,2770,2790,2810,2830,2850,1350,1400,1450,1500,700,800,1000,1100,1350,1400,900,1100,200,800,1066.67,1333.33,1600],'cm^-1')),
        HinderedRotor(inertia=(0.156089,'amu*angstrom^2'), symmetry=1, barrier=(3.5888,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.156089,'amu*angstrom^2'), symmetry=1, barrier=(3.5888,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.156089,'amu*angstrom^2'), symmetry=1, barrier=(3.5888,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.156089,'amu*angstrom^2'), symmetry=1, barrier=(3.5888,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.156089,'amu*angstrom^2'), symmetry=1, barrier=(3.5888,'kJ/mol'), semiclassical=False),
    ],
    spinMultiplicity = 3,
    opticalIsomers = 1,
    molecularWeight = (84.1595,'amu'),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[0.808226,0.0604888,-2.37075e-05,-6.10858e-09,5.08143e-12,30672.6,26.2777], Tmin=(100,'K'), Tmax=(1098.6,'K')), NASAPolynomial(coeffs=[12.2373,0.0362549,-1.43486e-05,2.6121e-09,-1.79974e-13,27112.7,-34.7092], Tmin=(1098.6,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(254.004,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(411.566,'J/(mol*K)'), comment="""Thermo group additivity estimation: group(Cs-CsCsCsH) + group(Cs-CsCsHH) + longDistanceInteraction_noncyclic(CsCs-ST) + group(Cs-CsCsHH) + longDistanceInteraction_noncyclic(CsCs-ST) + group(Cs-CsHHH) + group(Cs-CsHHH) + group(Cs-CsHHH) + radical(CCJ2_triplet)"""),
)

species(
    label = 'C=C(C)C(C)C(170)',
    structure = SMILES('C=C(C)C(C)C'),
    E0 = (-83.1449,'kJ/mol'),
    modes = [
        HarmonicOscillator(frequencies=([1380,1390,370,380,2900,435,2950,3100,1380,975,1025,1650,350,440,435,1725,2750,2762.5,2775,2787.5,2800,2812.5,2825,2837.5,2850,1350,1380,1410,1440,1470,1500,700,750,800,1000,1050,1100,1350,1375,1400,900,1000,1100,292.543],'cm^-1')),
        HinderedRotor(inertia=(0.00205755,'amu*angstrom^2'), symmetry=1, barrier=(0.119627,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.192175,'amu*angstrom^2'), symmetry=1, barrier=(11.9459,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.198364,'amu*angstrom^2'), symmetry=1, barrier=(11.9798,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.196892,'amu*angstrom^2'), symmetry=1, barrier=(11.9963,'kJ/mol'), semiclassical=False),
    ],
    spinMultiplicity = 1,
    opticalIsomers = 1,
    molecularWeight = (84.1595,'amu'),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[1.08853,0.0532258,-5.17688e-06,-2.43859e-08,1.1371e-11,-9885.93,22.593], Tmin=(100,'K'), Tmax=(1042.83,'K')), NASAPolynomial(coeffs=[11.4671,0.0364696,-1.42343e-05,2.60304e-09,-1.81103e-13,-13304,-33.9229], Tmin=(1042.83,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(-83.1449,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(415.724,'J/(mol*K)'), comment="""Thermo group additivity estimation: group(Cs-(Cds-Cds)CsCsH) + group(Cs-CsHHH) + group(Cs-CsHHH) + group(Cs-(Cds-Cds)HHH) + group(Cds-CdsCsCs) + longDistanceInteraction_noncyclic(CdCs-ST) + group(Cds-CdsHH)"""),
)

species(
    label = 'C[C](C)C(C)C(171)',
    structure = SMILES('C[C](C)C(C)C'),
    E0 = (-17.8305,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
    molecularWeight = (85.1674,'amu'),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[1.17493,0.0594329,-2.93088e-05,6.6791e-09,-5.94557e-13,-2040.42,23.452], Tmin=(100,'K'), Tmax=(2561.91,'K')), NASAPolynomial(coeffs=[22.9371,0.0254544,-9.41406e-06,1.50195e-09,-8.93443e-14,-13190.8,-102.01], Tmin=(2561.91,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(-17.8305,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(436.51,'J/(mol*K)'), comment="""Thermo group additivity estimation: group(Cs-CsCsCsH) + longDistanceInteraction_noncyclic(CsCs-TT) + group(Cs-CsCsCsH) + longDistanceInteraction_noncyclic(CsCs-TT) + group(Cs-CsHHH) + group(Cs-CsHHH) + group(Cs-CsHHH) + group(Cs-CsHHH) + radical(Tertalkyl)"""),
)

species(
    label = '[CH2]C(C)[C](C)C(172)',
    structure = SMILES('[CH2]C(C)[C](C)C'),
    E0 = (187.252,'kJ/mol'),
    modes = [
        HarmonicOscillator(frequencies=([1380,1390,370,380,2900,435,3000,3100,440,815,1455,1000,360,370,350,2750,2762.5,2775,2787.5,2800,2812.5,2825,2837.5,2850,1350,1380,1410,1440,1470,1500,700,750,800,1000,1050,1100,1350,1375,1400,900,1000,1100,180],'cm^-1')),
        HinderedRotor(inertia=(0.130717,'amu*angstrom^2'), symmetry=1, barrier=(3.00544,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.130406,'amu*angstrom^2'), symmetry=1, barrier=(2.9983,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.130349,'amu*angstrom^2'), symmetry=1, barrier=(2.99698,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(4.09962,'amu*angstrom^2'), symmetry=1, barrier=(94.2582,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.000285111,'amu*angstrom^2'), symmetry=1, barrier=(3.00615,'kJ/mol'), semiclassical=False),
    ],
    spinMultiplicity = 3,
    opticalIsomers = 1,
    molecularWeight = (84.1595,'amu'),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[1.48946,0.0581244,-3.34564e-05,1.09726e-08,-1.68144e-12,22608.8,24.7615], Tmin=(100,'K'), Tmax=(1311.91,'K')), NASAPolynomial(coeffs=[5.80137,0.0449775,-1.84246e-05,3.33395e-09,-2.2582e-13,21477.4,2.78838], Tmin=(1311.91,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(187.252,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(411.566,'J/(mol*K)'), comment="""Thermo group additivity estimation: group(Cs-CsCsCsH) + longDistanceInteraction_noncyclic(CsCs-TT) + group(Cs-CsCsCsH) + longDistanceInteraction_noncyclic(CsCs-TT) + group(Cs-CsHHH) + group(Cs-CsHHH) + group(Cs-CsHHH) + group(Cs-CsHHH) + radical(Isobutyl) + radical(Tertalkyl)"""),
)

species(
    label = '[CH2][C](C)C(C)C(173)',
    structure = SMILES('[CH2][C](C)C(C)C'),
    E0 = (187.252,'kJ/mol'),
    modes = [
        HarmonicOscillator(frequencies=([1380,1390,370,380,2900,435,3000,3100,440,815,1455,1000,360,370,350,2750,2762.5,2775,2787.5,2800,2812.5,2825,2837.5,2850,1350,1380,1410,1440,1470,1500,700,750,800,1000,1050,1100,1350,1375,1400,900,1000,1100,180],'cm^-1')),
        HinderedRotor(inertia=(0.130717,'amu*angstrom^2'), symmetry=1, barrier=(3.00544,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.130406,'amu*angstrom^2'), symmetry=1, barrier=(2.9983,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.130349,'amu*angstrom^2'), symmetry=1, barrier=(2.99698,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(4.09962,'amu*angstrom^2'), symmetry=1, barrier=(94.2582,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.000285111,'amu*angstrom^2'), symmetry=1, barrier=(3.00615,'kJ/mol'), semiclassical=False),
    ],
    spinMultiplicity = 3,
    opticalIsomers = 1,
    molecularWeight = (84.1595,'amu'),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[1.48946,0.0581244,-3.34564e-05,1.09726e-08,-1.68144e-12,22608.8,24.7615], Tmin=(100,'K'), Tmax=(1311.91,'K')), NASAPolynomial(coeffs=[5.80137,0.0449775,-1.84246e-05,3.33395e-09,-2.2582e-13,21477.4,2.78838], Tmin=(1311.91,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(187.252,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(411.566,'J/(mol*K)'), comment="""Thermo group additivity estimation: group(Cs-CsCsCsH) + longDistanceInteraction_noncyclic(CsCs-TT) + group(Cs-CsCsCsH) + longDistanceInteraction_noncyclic(CsCs-TT) + group(Cs-CsHHH) + group(Cs-CsHHH) + group(Cs-CsHHH) + group(Cs-CsHHH) + radical(Isobutyl) + radical(Tertalkyl)"""),
)

species(
    label = '[CH2]C(C)C([CH2])C(174)',
    structure = SMILES('[CH2]C(C)C([CH2])C'),
    E0 = (206.912,'kJ/mol'),
    modes = [
        HarmonicOscillator(frequencies=([1380,1383.33,1386.67,1390,370,373.333,376.667,380,2800,3000,430,440,3000,3033.33,3066.67,3100,415,465,780,850,1435,1475,900,1100,2750,2770,2790,2810,2830,2850,1350,1400,1450,1500,700,800,1000,1100,1350,1400,900,1100,1661.63],'cm^-1')),
        HinderedRotor(inertia=(0.00401062,'amu*angstrom^2'), symmetry=1, barrier=(0.119627,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.00520297,'amu*angstrom^2'), symmetry=1, barrier=(0.119627,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.407303,'amu*angstrom^2'), symmetry=1, barrier=(9.3647,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.00378038,'amu*angstrom^2'), symmetry=1, barrier=(0.119627,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(2.58568,'amu*angstrom^2'), symmetry=1, barrier=(69.9337,'kJ/mol'), semiclassical=False),
    ],
    spinMultiplicity = 3,
    opticalIsomers = 1,
    molecularWeight = (84.1595,'amu'),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[1.02611,0.0569757,-1.48991e-05,-1.81004e-08,1.09672e-11,25000.4,27.0772], Tmin=(100,'K'), Tmax=(951.683,'K')), NASAPolynomial(coeffs=[10.9871,0.0351464,-1.2075e-05,2.04503e-09,-1.36546e-13,22197,-25.2528], Tmin=(951.683,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(206.912,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(411.566,'J/(mol*K)'), comment="""Thermo group additivity estimation: group(Cs-CsCsCsH) + longDistanceInteraction_noncyclic(CsCs-TT) + group(Cs-CsCsCsH) + longDistanceInteraction_noncyclic(CsCs-TT) + group(Cs-CsHHH) + group(Cs-CsHHH) + group(Cs-CsHHH) + group(Cs-CsHHH) + radical(Isobutyl) + radical(Isobutyl)"""),
)

species(
    label = '[CH2]C([CH2])C(C)C(175)',
    structure = SMILES('[CH2]C([CH2])C(C)C'),
    E0 = (206.912,'kJ/mol'),
    modes = [
        HarmonicOscillator(frequencies=([1380,1383.33,1386.67,1390,370,373.333,376.667,380,2800,3000,430,440,3000,3033.33,3066.67,3100,415,465,780,850,1435,1475,900,1100,2750,2770,2790,2810,2830,2850,1350,1400,1450,1500,700,800,1000,1100,1350,1400,900,1100,1661.63],'cm^-1')),
        HinderedRotor(inertia=(0.00401062,'amu*angstrom^2'), symmetry=1, barrier=(0.119627,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.00520297,'amu*angstrom^2'), symmetry=1, barrier=(0.119627,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.407303,'amu*angstrom^2'), symmetry=1, barrier=(9.3647,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.00378038,'amu*angstrom^2'), symmetry=1, barrier=(0.119627,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(2.58568,'amu*angstrom^2'), symmetry=1, barrier=(69.9337,'kJ/mol'), semiclassical=False),
    ],
    spinMultiplicity = 3,
    opticalIsomers = 1,
    molecularWeight = (84.1595,'amu'),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[1.02611,0.0569757,-1.48991e-05,-1.81004e-08,1.09672e-11,25000.4,26.3841], Tmin=(100,'K'), Tmax=(951.683,'K')), NASAPolynomial(coeffs=[10.9871,0.0351464,-1.2075e-05,2.04503e-09,-1.36546e-13,22197,-25.9459], Tmin=(951.683,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(206.912,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(411.566,'J/(mol*K)'), comment="""Thermo group additivity estimation: group(Cs-CsCsCsH) + longDistanceInteraction_noncyclic(CsCs-TT) + group(Cs-CsCsCsH) + longDistanceInteraction_noncyclic(CsCs-TT) + group(Cs-CsHHH) + group(Cs-CsHHH) + group(Cs-CsHHH) + group(Cs-CsHHH) + radical(Isobutyl) + radical(Isobutyl)"""),
)

species(
    label = '[CH]C(C)C(C)C(176)',
    structure = SMILES('[CH]C(C)C(C)C'),
    E0 = (244.962,'kJ/mol'),
    modes = [
        HarmonicOscillator(frequencies=([1380,1383.33,1386.67,1390,370,373.333,376.667,380,2800,3000,430,440,2750,2762.5,2775,2787.5,2800,2812.5,2825,2837.5,2850,1350,1380,1410,1440,1470,1500,700,750,800,1000,1050,1100,1350,1375,1400,900,1000,1100,200,800,1200,1600],'cm^-1')),
        HinderedRotor(inertia=(0.156089,'amu*angstrom^2'), symmetry=1, barrier=(3.5888,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.156089,'amu*angstrom^2'), symmetry=1, barrier=(3.5888,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.156089,'amu*angstrom^2'), symmetry=1, barrier=(3.5888,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.156089,'amu*angstrom^2'), symmetry=1, barrier=(3.5888,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.156089,'amu*angstrom^2'), symmetry=1, barrier=(3.5888,'kJ/mol'), semiclassical=False),
    ],
    spinMultiplicity = 3,
    opticalIsomers = 1,
    molecularWeight = (84.1595,'amu'),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[0.892792,0.0571654,-1.10585e-05,-2.14445e-08,1.10298e-11,29583.7,24.998], Tmin=(100,'K'), Tmax=(1031.03,'K')), NASAPolynomial(coeffs=[12.6465,0.0352434,-1.36132e-05,2.48171e-09,-1.72743e-13,25901.5,-38.169], Tmin=(1031.03,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(244.962,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(411.566,'J/(mol*K)'), comment="""Thermo group additivity estimation: group(Cs-CsCsCsH) + longDistanceInteraction_noncyclic(CsCs-TT) + group(Cs-CsCsCsH) + longDistanceInteraction_noncyclic(CsCs-TT) + group(Cs-CsHHH) + group(Cs-CsHHH) + group(Cs-CsHHH) + group(Cs-CsHHH) + radical(CCJ2_triplet)"""),
)

species(
    label = 'C=CC(C)(C)C(177)',
    structure = SMILES('C=CC(C)(C)C'),
    E0 = (-80.9798,'kJ/mol'),
    modes = [
        HarmonicOscillator(frequencies=([350,400,1190,1240,450,2950,3100,1380,975,1025,1650,3010,987.5,1337.5,450,1655,2750,2762.5,2775,2787.5,2800,2812.5,2825,2837.5,2850,1350,1380,1410,1440,1470,1500,700,750,800,1000,1050,1100,1350,1375,1400,900,1000,1100,180],'cm^-1')),
        HinderedRotor(inertia=(0.6948,'amu*angstrom^2'), symmetry=1, barrier=(15.9748,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.693572,'amu*angstrom^2'), symmetry=1, barrier=(15.9466,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.694733,'amu*angstrom^2'), symmetry=1, barrier=(15.9733,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.693798,'amu*angstrom^2'), symmetry=1, barrier=(15.9518,'kJ/mol'), semiclassical=False),
    ],
    spinMultiplicity = 1,
    opticalIsomers = 1,
    molecularWeight = (84.1595,'amu'),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[1.11139,0.0483501,1.78644e-05,-5.2716e-08,2.20236e-11,-9622.33,20.8354], Tmin=(100,'K'), Tmax=(1012.16,'K')), NASAPolynomial(coeffs=[13.7443,0.0339319,-1.33873e-05,2.526e-09,-1.81283e-13,-13998.4,-49.2482], Tmin=(1012.16,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(-80.9798,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(415.724,'J/(mol*K)'), comment="""Thermo group additivity estimation: group(Cs-(Cds-Cds)CsCsCs) + group(Cs-CsHHH) + group(Cs-CsHHH) + group(Cs-CsHHH) + group(Cds-CdsCsH) + group(Cds-CdsHH)"""),
)

species(
    label = 'CC=C(C)C(178)',
    structure = SMILES('CC=C(C)C'),
    E0 = (-68.9569,'kJ/mol'),
    modes = [
        HarmonicOscillator(frequencies=([2750,2762.5,2775,2787.5,2800,2812.5,2825,2837.5,2850,1350,1380,1410,1440,1470,1500,700,750,800,1000,1050,1100,1350,1375,1400,900,1000,1100,3010,987.5,1337.5,450,1655,350,440,435,1725],'cm^-1')),
        HinderedRotor(inertia=(0.415137,'amu*angstrom^2'), symmetry=1, barrier=(9.54483,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.414794,'amu*angstrom^2'), symmetry=1, barrier=(9.53693,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.414559,'amu*angstrom^2'), symmetry=1, barrier=(9.53153,'kJ/mol'), semiclassical=False),
    ],
    spinMultiplicity = 1,
    opticalIsomers = 1,
    molecularWeight = (70.1329,'amu'),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[1.7164,0.0447584,-1.90576e-05,1.71722e-09,5.44318e-13,-8206.95,17.9461], Tmin=(100,'K'), Tmax=(1463.13,'K')), NASAPolynomial(coeffs=[9.80821,0.0303176,-1.21276e-05,2.14766e-09,-1.42315e-13,-11397,-26.9814], Tmin=(1463.13,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(-68.9569,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(345.051,'J/(mol*K)'), comment="""Thermo group additivity estimation: group(Cs-(Cds-Cds)HHH) + group(Cs-(Cds-Cds)HHH) + group(Cs-(Cds-Cds)HHH) + group(Cds-CdsCsCs) + group(Cds-CdsCsH)"""),
)

species(
    label = 'C[CH][C](C)C(179)',
    structure = SMILES('C[CH][C](C)C'),
    E0 = (206.186,'kJ/mol'),
    modes = [
        HarmonicOscillator(frequencies=([2750,2762.5,2775,2787.5,2800,2812.5,2825,2837.5,2850,1350,1380,1410,1440,1470,1500,700,750,800,1000,1050,1100,1350,1375,1400,900,1000,1100,3025,407.5,1350,352.5,360,370,350,2286.11],'cm^-1')),
        HinderedRotor(inertia=(0.0737563,'amu*angstrom^2'), symmetry=1, barrier=(1.69749,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.302096,'amu*angstrom^2'), symmetry=1, barrier=(6.94649,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.00187647,'amu*angstrom^2'), symmetry=1, barrier=(6.9485,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.273525,'amu*angstrom^2'), symmetry=1, barrier=(6.94587,'kJ/mol'), semiclassical=False),
    ],
    spinMultiplicity = 3,
    opticalIsomers = 1,
    molecularWeight = (70.1329,'amu'),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[2.67808,0.0384653,-1.27677e-05,1.43436e-10,3.68101e-13,24837.4,19.3023], Tmin=(100,'K'), Tmax=(2101.54,'K')), NASAPolynomial(coeffs=[15.6717,0.0243403,-1.02564e-05,1.74837e-09,-1.08515e-13,17033.8,-58.6069], Tmin=(2101.54,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(206.186,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(340.893,'J/(mol*K)'), comment="""Thermo group additivity estimation: group(Cs-CsCsCsH) + group(Cs-CsCsHH) + longDistanceInteraction_noncyclic(CsCs-ST) + group(Cs-CsHHH) + group(Cs-CsHHH) + group(Cs-CsHHH) + radical(Cs_S) + radical(Tertalkyl)"""),
)

species(
    label = '[CH2][CH]C(C)(C)C(180)',
    structure = SMILES('[CH2][CH]C(C)(C)C'),
    E0 = (193.517,'kJ/mol'),
    modes = [
        HarmonicOscillator(frequencies=([350,400,1190,1240,450,3025,407.5,1350,352.5,3000,3100,440,815,1455,1000,2750,2762.5,2775,2787.5,2800,2812.5,2825,2837.5,2850,1350,1380,1410,1440,1470,1500,700,750,800,1000,1050,1100,1350,1375,1400,900,1000,1100,851.977],'cm^-1')),
        HinderedRotor(inertia=(0.113915,'amu*angstrom^2'), symmetry=1, barrier=(2.61914,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.114029,'amu*angstrom^2'), symmetry=1, barrier=(2.62175,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.113993,'amu*angstrom^2'), symmetry=1, barrier=(2.62092,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.657472,'amu*angstrom^2'), symmetry=1, barrier=(15.1166,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.657425,'amu*angstrom^2'), symmetry=1, barrier=(15.1155,'kJ/mol'), semiclassical=False),
    ],
    spinMultiplicity = 3,
    opticalIsomers = 1,
    molecularWeight = (84.1595,'amu'),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[0.984454,0.0558635,-1.33135e-05,-1.43027e-08,7.1044e-12,23391.9,25.1965], Tmin=(100,'K'), Tmax=(1132.62,'K')), NASAPolynomial(coeffs=[12.0445,0.0371869,-1.55744e-05,2.91774e-09,-2.03916e-13,19579.1,-35.3111], Tmin=(1132.62,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(193.517,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(411.566,'J/(mol*K)'), comment="""Thermo group additivity estimation: group(Cs-CsCsCsCs) + group(Cs-CsCsHH) + longDistanceInteraction_noncyclic(CsCs-SQ) + group(Cs-CsHHH) + group(Cs-CsHHH) + group(Cs-CsHHH) + group(Cs-CsHHH) + radical(RCCJ) + radical(Cs_S)"""),
)

species(
    label = 'C[C](C)C(181)',
    structure = SMILES('C[C](C)C'),
    E0 = (32.0771,'kJ/mol'),
    modes = [
        HarmonicOscillator(frequencies=([2750,2762.5,2775,2787.5,2800,2812.5,2825,2837.5,2850,1350,1380,1410,1440,1470,1500,700,750,800,1000,1050,1100,1350,1375,1400,900,1000,1100,360,370,350],'cm^-1')),
        HinderedRotor(inertia=(0.0936296,'amu*angstrom^2'), symmetry=1, barrier=(2.15273,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.0936438,'amu*angstrom^2'), symmetry=1, barrier=(2.15305,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.0938508,'amu*angstrom^2'), symmetry=1, barrier=(2.15782,'kJ/mol'), semiclassical=False),
    ],
    spinMultiplicity = 2,
    opticalIsomers = 1,
    molecularWeight = (57.1143,'amu'),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[2.58814,0.0309399,-7.08916e-06,-2.52748e-09,9.71938e-13,3908.8,13.7543], Tmin=(100,'K'), Tmax=(1673.54,'K')), NASAPolynomial(coeffs=[7.58069,0.0271246,-1.09452e-05,1.90694e-09,-1.23457e-13,1100.99,-16.2991], Tmin=(1673.54,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(32.0771,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(295.164,'J/(mol*K)'), comment="""Thermo group additivity estimation: group(Cs-CsCsCsH) + group(Cs-CsHHH) + group(Cs-CsHHH) + group(Cs-CsHHH) + radical(Tertalkyl)"""),
)

species(
    label = '[CH]C(C)(C)C(182)',
    structure = SMILES('[CH]C(C)(C)C'),
    E0 = (259.029,'kJ/mol'),
    modes = [
        HarmonicOscillator(frequencies=([350,400,1190,1240,450,2750,2762.5,2775,2787.5,2800,2812.5,2825,2837.5,2850,1350,1380,1410,1440,1470,1500,700,750,800,1000,1050,1100,1350,1375,1400,900,1000,1100,180,426.391,4000],'cm^-1')),
        HinderedRotor(inertia=(0.836042,'amu*angstrom^2'), symmetry=1, barrier=(19.2222,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.0155331,'amu*angstrom^2'), symmetry=1, barrier=(19.2406,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.835792,'amu*angstrom^2'), symmetry=1, barrier=(19.2165,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.834958,'amu*angstrom^2'), symmetry=1, barrier=(19.1973,'kJ/mol'), semiclassical=False),
    ],
    spinMultiplicity = 3,
    opticalIsomers = 1,
    molecularWeight = (70.1329,'amu'),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[1.47517,0.0432754,6.93669e-06,-3.74329e-08,1.65795e-11,31255.7,17.8434], Tmin=(100,'K'), Tmax=(1006.46,'K')), NASAPolynomial(coeffs=[12.8021,0.026557,-1.03218e-05,1.93517e-09,-1.38643e-13,27542.4,-43.9956], Tmin=(1006.46,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(259.029,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(340.893,'J/(mol*K)'), comment="""Thermo group additivity estimation: group(Cs-CsCsCsCs) + group(Cs-CsHHH) + group(Cs-CsHHH) + group(Cs-CsHHH) + group(Cs-CsHHH) + radical(CCJ2_triplet)"""),
)

species(
    label = 'C[C]C(C)(C)C(183)',
    structure = SMILES('C[C]C(C)(C)C'),
    E0 = (241.944,'kJ/mol'),
    modes = [
        HarmonicOscillator(frequencies=([350,400,1190,1240,450,2750,2759.09,2768.18,2777.27,2786.36,2795.45,2804.55,2813.64,2822.73,2831.82,2840.91,2850,1350,1371.43,1392.86,1414.29,1435.71,1457.14,1478.57,1500,700,733.333,766.667,800,1000,1033.33,1066.67,1100,1350,1366.67,1383.33,1400,900,966.667,1033.33,1100,180,662.938],'cm^-1')),
        HinderedRotor(inertia=(0.695935,'amu*angstrom^2'), symmetry=1, barrier=(16.0009,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.697451,'amu*angstrom^2'), symmetry=1, barrier=(16.0358,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.695632,'amu*angstrom^2'), symmetry=1, barrier=(15.994,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.696662,'amu*angstrom^2'), symmetry=1, barrier=(16.0176,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.696282,'amu*angstrom^2'), symmetry=1, barrier=(16.0089,'kJ/mol'), semiclassical=False),
    ],
    spinMultiplicity = 3,
    opticalIsomers = 1,
    molecularWeight = (84.1595,'amu'),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[0.728645,0.0592265,-1.09334e-05,-2.42486e-08,1.23035e-11,29227.8,21.6703], Tmin=(100,'K'), Tmax=(1039.07,'K')), NASAPolynomial(coeffs=[14.341,0.033928,-1.35393e-05,2.52714e-09,-1.78676e-13,24935.9,-51.564], Tmin=(1039.07,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(241.944,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(411.566,'J/(mol*K)'), comment="""Thermo group additivity estimation: group(Cs-CsCsCsCs) + group(Cs-CsCsHH) + longDistanceInteraction_noncyclic(CsCs-SQ) + group(Cs-CsHHH) + group(Cs-CsHHH) + group(Cs-CsHHH) + group(Cs-CsHHH) + radical(CCJ2_triplet)"""),
)

species(
    label = '[CH]CC(C)(C)C(184)',
    structure = SMILES('[CH]CC(C)(C)C'),
    E0 = (241.944,'kJ/mol'),
    modes = [
        HarmonicOscillator(frequencies=([350,400,1190,1240,450,2750,2850,1437.5,1250,1305,750,350,2750,2762.5,2775,2787.5,2800,2812.5,2825,2837.5,2850,1350,1380,1410,1440,1470,1500,700,750,800,1000,1050,1100,1350,1375,1400,900,1000,1100,200,800,1200,1600],'cm^-1')),
        HinderedRotor(inertia=(0.156089,'amu*angstrom^2'), symmetry=1, barrier=(3.5888,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.156089,'amu*angstrom^2'), symmetry=1, barrier=(3.5888,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.156089,'amu*angstrom^2'), symmetry=1, barrier=(3.5888,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.156089,'amu*angstrom^2'), symmetry=1, barrier=(3.5888,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.156089,'amu*angstrom^2'), symmetry=1, barrier=(3.5888,'kJ/mol'), semiclassical=False),
    ],
    spinMultiplicity = 3,
    opticalIsomers = 1,
    molecularWeight = (84.1595,'amu'),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[0.728645,0.0592265,-1.09334e-05,-2.42486e-08,1.23035e-11,29227.8,22.7689], Tmin=(100,'K'), Tmax=(1039.07,'K')), NASAPolynomial(coeffs=[14.341,0.033928,-1.35393e-05,2.52714e-09,-1.78676e-13,24935.9,-50.4653], Tmin=(1039.07,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(241.944,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(411.566,'J/(mol*K)'), comment="""Thermo group additivity estimation: group(Cs-CsCsCsCs) + group(Cs-CsCsHH) + longDistanceInteraction_noncyclic(CsCs-SQ) + group(Cs-CsHHH) + group(Cs-CsHHH) + group(Cs-CsHHH) + group(Cs-CsHHH) + radical(CCJ2_triplet)"""),
)

species(
    label = 'Ne',
    structure = SMILES('[Ne]'),
    E0 = (-6.19738,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
    molecularWeight = (20.1797,'amu'),
    collisionModel = TransportData(shapeIndex=0, epsilon=(1235.53,'J/mol'), sigma=(3.758e-10,'m'), dipoleMoment=(0,'C*m'), polarizability=(0,'angstroms^3'), rotrelaxcollnum=0, comment="""Epsilon & sigma estimated with fixed Lennard Jones Parameters. This is the fallback method! Try improving transport databases!"""),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[2.5,0,0,0,0,-745.375,3.35532], Tmin=(200,'K'), Tmax=(1000,'K')), NASAPolynomial(coeffs=[2.5,0,0,0,0,-745.375,3.35532], Tmin=(1000,'K'), Tmax=(6000,'K'))], Tmin=(200,'K'), Tmax=(6000,'K'), E0=(-6.19738,'kJ/mol'), Cp0=(20.7862,'J/(mol*K)'), CpInf=(20.7862,'J/(mol*K)'), label="""Ne""", comment="""Thermo library: primaryThermoLibrary"""),
)

species(
    label = 'Ar',
    structure = SMILES('[Ar]'),
    E0 = (-6.19738,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
    molecularWeight = (39.348,'amu'),
    collisionModel = TransportData(shapeIndex=0, epsilon=(1134.93,'J/mol'), sigma=(3.33,'angstroms'), dipoleMoment=(0,'C*m'), polarizability=(0,'angstroms^3'), rotrelaxcollnum=0.0, comment="""GRI-Mech"""),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[2.5,0,0,0,0,-745.375,4.37967], Tmin=(200,'K'), Tmax=(1000,'K')), NASAPolynomial(coeffs=[2.5,0,0,0,0,-745.375,4.37967], Tmin=(1000,'K'), Tmax=(6000,'K'))], Tmin=(200,'K'), Tmax=(6000,'K'), E0=(-6.19738,'kJ/mol'), Cp0=(20.7862,'J/(mol*K)'), CpInf=(20.7862,'J/(mol*K)'), label="""Ar""", comment="""Thermo library: primaryThermoLibrary"""),
)

species(
    label = 'He',
    structure = SMILES('[He]'),
    E0 = (-6.19738,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
    molecularWeight = (4.0026,'amu'),
    collisionModel = TransportData(shapeIndex=0, epsilon=(84.8076,'J/mol'), sigma=(2.576,'angstroms'), dipoleMoment=(0,'De'), polarizability=(0,'angstroms^3'), rotrelaxcollnum=0.0, comment="""NOx2018"""),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[2.5,0,0,0,0,-745.375,0.928724], Tmin=(200,'K'), Tmax=(1000,'K')), NASAPolynomial(coeffs=[2.5,0,0,0,0,-745.375,0.928724], Tmin=(1000,'K'), Tmax=(6000,'K'))], Tmin=(200,'K'), Tmax=(6000,'K'), E0=(-6.19738,'kJ/mol'), Cp0=(20.7862,'J/(mol*K)'), CpInf=(20.7862,'J/(mol*K)'), label="""He""", comment="""Thermo library: primaryThermoLibrary"""),
)

species(
    label = 'N2',
    structure = SMILES('N#N'),
    E0 = (-8.64289,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
    molecularWeight = (28.0135,'amu'),
    collisionModel = TransportData(shapeIndex=1, epsilon=(810.913,'J/mol'), sigma=(3.621,'angstroms'), dipoleMoment=(0,'C*m'), polarizability=(1.76,'angstroms^3'), rotrelaxcollnum=4.0, comment="""PrimaryTransportLibrary"""),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[3.53101,-0.000123661,-5.02999e-07,2.43531e-09,-1.40881e-12,-1046.98,2.96747], Tmin=(200,'K'), Tmax=(1000,'K')), NASAPolynomial(coeffs=[2.95258,0.0013969,-4.92632e-07,7.8601e-11,-4.60755e-15,-923.949,5.87189], Tmin=(1000,'K'), Tmax=(6000,'K'))], Tmin=(200,'K'), Tmax=(6000,'K'), E0=(-8.64289,'kJ/mol'), Cp0=(29.1007,'J/(mol*K)'), CpInf=(37.4151,'J/(mol*K)'), label="""N2""", comment="""Thermo library: primaryThermoLibrary"""),
)

transitionState(
    label = 'TS1',
    E0 = (363.915,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS2',
    E0 = (141.209,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS3',
    E0 = (150.233,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS4',
    E0 = (110.72,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS5',
    E0 = (161.586,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS6',
    E0 = (59.8242,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS7',
    E0 = (364.171,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS8',
    E0 = (408.845,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS9',
    E0 = (367.076,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS10',
    E0 = (408.833,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS11',
    E0 = (408.845,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS12',
    E0 = (419.633,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS13',
    E0 = (419.633,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS14',
    E0 = (445.453,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS15',
    E0 = (167.459,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS16',
    E0 = (446.73,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS17',
    E0 = (449.828,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS18',
    E0 = (468.155,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS19',
    E0 = (157.764,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS20',
    E0 = (124.263,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS21',
    E0 = (154.801,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS22',
    E0 = (374.692,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS23',
    E0 = (374.715,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS24',
    E0 = (419.645,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS25',
    E0 = (375.141,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS26',
    E0 = (419.645,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS27',
    E0 = (377.876,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS28',
    E0 = (419.633,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS29',
    E0 = (419.633,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS30',
    E0 = (430.433,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS31',
    E0 = (456.253,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS32',
    E0 = (452.848,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS33',
    E0 = (468.155,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS34',
    E0 = (142.344,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS35',
    E0 = (141.209,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS36',
    E0 = (111.879,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS37',
    E0 = (118.154,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS38',
    E0 = (364.318,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS39',
    E0 = (367.088,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS40',
    E0 = (408.845,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS41',
    E0 = (408.857,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS42',
    E0 = (366.99,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS43',
    E0 = (408.845,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS44',
    E0 = (419.645,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS45',
    E0 = (419.645,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS46',
    E0 = (445.465,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS47',
    E0 = (445.453,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS48',
    E0 = (170.806,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS49',
    E0 = (202.289,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS50',
    E0 = (446.73,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS51',
    E0 = (446.753,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS52',
    E0 = (468.155,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS53',
    E0 = (132.448,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS54',
    E0 = (118.21,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS55',
    E0 = (124.828,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS56',
    E0 = (110.203,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS57',
    E0 = (148.943,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS58',
    E0 = (94.1323,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS59',
    E0 = (60.4059,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS60',
    E0 = (364.171,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS61',
    E0 = (366.99,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS62',
    E0 = (404.751,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS63',
    E0 = (365.346,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS64',
    E0 = (413.87,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS65',
    E0 = (372.017,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS66',
    E0 = (413.774,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS67',
    E0 = (424.411,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS68',
    E0 = (424.575,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS69',
    E0 = (456.253,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS70',
    E0 = (450.395,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS71',
    E0 = (442.048,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS72',
    E0 = (462.461,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS73',
    E0 = (123.759,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS74',
    E0 = (132.448,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS75',
    E0 = (93.8373,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS76',
    E0 = (151.373,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS77',
    E0 = (141.653,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS78',
    E0 = (117.267,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS79',
    E0 = (346.089,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS80',
    E0 = (352.358,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS81',
    E0 = (394.115,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS82',
    E0 = (394.211,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS83',
    E0 = (404.915,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS84',
    E0 = (404.751,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS85',
    E0 = (430.735,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS86',
    E0 = (156.045,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS87',
    E0 = (446.73,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS88',
    E0 = (449.828,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS89',
    E0 = (148.723,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS90',
    E0 = (115.057,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS91',
    E0 = (155.383,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS92',
    E0 = (158.485,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS93',
    E0 = (363.769,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS94',
    E0 = (366.99,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS95',
    E0 = (404.915,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS96',
    E0 = (365.936,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS97',
    E0 = (414.034,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS98',
    E0 = (413.939,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS99',
    E0 = (424.575,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS100',
    E0 = (456.253,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS101',
    E0 = (447.154,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS102',
    E0 = (462.461,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS103',
    E0 = (129.282,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS104',
    E0 = (133.237,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS105',
    E0 = (111.285,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS106',
    E0 = (110.182,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS107',
    E0 = (155.975,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS108',
    E0 = (356.202,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS109',
    E0 = (394.211,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS110',
    E0 = (361.391,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS111',
    E0 = (403.234,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS112',
    E0 = (413.87,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS113',
    E0 = (414.034,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS114',
    E0 = (445.465,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS115',
    E0 = (439.855,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS116',
    E0 = (195.096,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS117',
    E0 = (196.594,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS118',
    E0 = (435.807,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS119',
    E0 = (437.712,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS120',
    E0 = (462.461,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS121',
    E0 = (133.582,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS122',
    E0 = (119.368,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS123',
    E0 = (113.55,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS124',
    E0 = (157.729,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS125',
    E0 = (93.7139,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS126',
    E0 = (364.318,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS127',
    E0 = (408.098,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS128',
    E0 = (371.853,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS129',
    E0 = (417.218,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS130',
    E0 = (427.922,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS131',
    E0 = (450.395,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS132',
    E0 = (442.06,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS133',
    E0 = (465.808,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS134',
    E0 = (133.237,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS135',
    E0 = (141.191,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS136',
    E0 = (94.7768,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS137',
    E0 = (352.846,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS138',
    E0 = (356.19,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS139',
    E0 = (394.115,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS140',
    E0 = (403.234,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS141',
    E0 = (413.774,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS142',
    E0 = (413.939,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS143',
    E0 = (445.453,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS144',
    E0 = (161.764,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS145',
    E0 = (437.524,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS146',
    E0 = (444.134,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS147',
    E0 = (462.461,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS148',
    E0 = (123.759,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS149',
    E0 = (133.582,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS150',
    E0 = (100.432,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS151',
    E0 = (154.721,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS152',
    E0 = (148.689,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS153',
    E0 = (352.108,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS154',
    E0 = (397.558,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS155',
    E0 = (408.262,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS156',
    E0 = (408.098,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS157',
    E0 = (430.735,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS158',
    E0 = (156.045,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS159',
    E0 = (446.753,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS160',
    E0 = (449.828,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS161',
    E0 = (129.282,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS162',
    E0 = (142.261,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS163',
    E0 = (104.691,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS164',
    E0 = (111.285,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS165',
    E0 = (117.579,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS166',
    E0 = (353.518,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS167',
    E0 = (356.202,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS168',
    E0 = (397.558,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS169',
    E0 = (361.313,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS170',
    E0 = (406.677,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS171',
    E0 = (417.218,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS172',
    E0 = (417.382,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS173',
    E0 = (417.382,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS174',
    E0 = (445.453,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS175',
    E0 = (439.855,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS176',
    E0 = (435.93,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS177',
    E0 = (444.134,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS178',
    E0 = (465.808,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS179',
    E0 = (106.85,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS180',
    E0 = (113.444,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS181',
    E0 = (145.585,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS182',
    E0 = (81.818,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS183',
    E0 = (346.089,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS184',
    E0 = (352.108,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS185',
    E0 = (356.278,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS186',
    E0 = (405.073,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS187',
    E0 = (415.529,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS188',
    E0 = (415.778,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS189',
    E0 = (450.395,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS190',
    E0 = (434.903,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS191',
    E0 = (427.33,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS192',
    E0 = (453.748,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS193',
    E0 = (149.792,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS194',
    E0 = (108.484,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS195',
    E0 = (363.915,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS196',
    E0 = (363.939,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS197',
    E0 = (367.002,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS198',
    E0 = (408.262,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS199',
    E0 = (372.017,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS200',
    E0 = (417.382,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS201',
    E0 = (417.382,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS202',
    E0 = (427.922,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS203',
    E0 = (428.086,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS204',
    E0 = (456.253,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS205',
    E0 = (450.559,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS206',
    E0 = (446.989,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS207',
    E0 = (465.808,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS208',
    E0 = (128.66,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS209',
    E0 = (102.308,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS210',
    E0 = (116.856,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS211',
    E0 = (104.509,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS212',
    E0 = (352.846,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS213',
    E0 = (361.227,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS214',
    E0 = (399.057,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS215',
    E0 = (361.391,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS216',
    E0 = (399.057,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS217',
    E0 = (418.716,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS218',
    E0 = (418.716,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS219',
    E0 = (450.395,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS220',
    E0 = (450.559,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS221',
    E0 = (436.449,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS222',
    E0 = (456.767,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS223',
    E0 = (133.335,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS224',
    E0 = (98.7669,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS225',
    E0 = (105.519,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS226',
    E0 = (341.568,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS227',
    E0 = (405.073,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS228',
    E0 = (405.322,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS229',
    E0 = (439.855,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS230',
    E0 = (145.589,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS231',
    E0 = (417.864,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS232',
    E0 = (428.727,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS233',
    E0 = (453.748,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS234',
    E0 = (140.867,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS235',
    E0 = (84.4773,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS236',
    E0 = (345.874,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS237',
    E0 = (352.272,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS238',
    E0 = (405.322,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS239',
    E0 = (415.778,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS240',
    E0 = (450.559,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS241',
    E0 = (431.498,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS242',
    E0 = (453.748,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

reaction(
    label = 'reaction21',
    reactants = ['[CH2]C[CH]C(44)', 'C[CH2](6)'],
    products = ['C[CH]CCCC(49)'],
    transitionState = 'TS1',
    kinetics = Arrhenius(A=(8.73e+14,'cm^3/(mol*s)'), n=-0.699, Ea=(0,'kJ/mol'), T0=(1,'K'), Tmin=(200,'K'), Tmax=(2000,'K'), comment="""From training reaction 11 used for C_rad/H2/Cs;C_rad/H2/Cs
Exact match found for rate rule [C_rad/H2/Cs;C_rad/H2/Cs]
Euclidian distance = 0
family: R_Recombination
Ea raised from -0.0 to 0 kJ/mol."""),
)

reaction(
    label = 'reaction17',
    reactants = ['H(8)', 'CC=CCCC(64)'],
    products = ['C[CH]CCCC(49)'],
    transitionState = 'TS2',
    kinetics = Arrhenius(A=(1.46e+08,'cm^3/(mol*s)'), n=1.64, Ea=(5.73208,'kJ/mol'), T0=(1,'K'), Tmin=(300,'K'), Tmax=(1500,'K'), comment="""From training reaction 2555 used for Cds-CsH_Cds-CsH;HJ
Exact match found for rate rule [Cds-CsH_Cds-CsH;HJ]
Euclidian distance = 0
family: R_Addition_MultipleBond"""),
)

reaction(
    label = 'reaction18',
    reactants = ['C=CCCCC(48)', 'H(8)'],
    products = ['C[CH]CCCC(49)'],
    transitionState = 'TS3',
    kinetics = Arrhenius(A=(3.36e+08,'cm^3/(mol*s)'), n=1.56, Ea=(2.5104,'kJ/mol'), T0=(1,'K'), comment="""From training reaction 9 used for Cds-HH_Cds-CsH;HJ
Exact match found for rate rule [Cds-HH_Cds-CsH;HJ]
Euclidian distance = 0
family: R_Addition_MultipleBond"""),
)

reaction(
    label = 'reaction19',
    reactants = ['[CH2]CC(4)', 'C=CC(36)'],
    products = ['C[CH]CCCC(49)'],
    transitionState = 'TS4',
    kinetics = Arrhenius(A=(2130,'cm^3/(mol*s)'), n=2.41, Ea=(19.874,'kJ/mol'), T0=(1,'K'), Tmin=(300,'K'), Tmax=(1500,'K'), comment="""From training reaction 238 used for Cds-HH_Cds-Cs\H3/H;CsJ-CsHH
Exact match found for rate rule [Cds-HH_Cds-Cs\H3/H;CsJ-CsHH]
Euclidian distance = 0
family: R_Addition_MultipleBond"""),
)

reaction(
    label = 'reaction5',
    reactants = ['CC[CH]CCC(50)'],
    products = ['C[CH]CCCC(49)'],
    transitionState = 'TS5',
    kinetics = Arrhenius(A=(6.76e+09,'s^-1'), n=0.88, Ea=(158.992,'kJ/mol'), T0=(1,'K'), Tmin=(300,'K'), Tmax=(1500,'K'), comment="""From training reaction 357 used for R2H_S;C_rad_out_H/NonDeC;Cs_H_out_H/NonDeC
Exact match found for rate rule [R2H_S;C_rad_out_H/NonDeC;Cs_H_out_H/NonDeC]
Euclidian distance = 0
Multiplied by reaction path degeneracy 2.0
family: intra_H_migration"""),
)

reaction(
    label = 'reaction6',
    reactants = ['[CH2]CCCCC(10)'],
    products = ['C[CH]CCCC(49)'],
    transitionState = 'TS6',
    kinetics = Arrhenius(A=(262000,'s^-1'), n=1.62, Ea=(11.1,'kcal/mol'), T0=(1,'K'), comment="""Matched reaction 112 CCCCC[CH2]-1 <=> C[CH]CCCC in intra_H_migration/training
This reaction matched rate rule [R2H_S;C_rad_out_H/NonDeC;Cs_H_out_2H]
family: intra_H_migration"""),
)

reaction(
    label = 'reaction22',
    reactants = ['[CH2]CC(4)', '[CH2][CH]C(38)'],
    products = ['C[CH]CCCC(49)'],
    transitionState = 'TS7',
    kinetics = Arrhenius(A=(6.25357e+06,'m^3/(mol*s)'), n=0.093384, Ea=(0.402701,'kJ/mol'), T0=(1,'K'), comment="""Estimated using an average for rate rule [Y_rad;C_rad/H2/Cs]
Euclidian distance = 0
family: R_Recombination"""),
)

reaction(
    label = 'reaction23',
    reactants = ['H(8)', 'C[CH]C[CH]CC(65)'],
    products = ['C[CH]CCCC(49)'],
    transitionState = 'TS8',
    kinetics = Arrhenius(A=(2e+13,'cm^3/(mol*s)','*|/',3.16), n=0, Ea=(0,'kJ/mol'), T0=(1,'K'), Tmin=(300,'K'), Tmax=(2000,'K'), comment="""From training reaction 59 used for H_rad;C_rad/H/NonDeC
Exact match found for rate rule [C_rad/H/NonDeC;H_rad]
Euclidian distance = 0
family: R_Recombination"""),
)

reaction(
    label = 'reaction24',
    reactants = ['[CH2]CC[CH]C(30)', '[CH3](11)'],
    products = ['C[CH]CCCC(49)'],
    transitionState = 'TS9',
    kinetics = Arrhenius(A=(1.23e+15,'cm^3/(mol*s)'), n=-0.562, Ea=(0.085772,'kJ/mol'), T0=(1,'K'), Tmin=(200,'K'), Tmax=(2000,'K'), comment="""From training reaction 10 used for C_methyl;C_rad/H2/Cs
Exact match found for rate rule [C_rad/H2/Cs;C_methyl]
Euclidian distance = 0
family: R_Recombination"""),
)

reaction(
    label = 'reaction25',
    reactants = ['C[CH]CC[CH]C(66)', 'H(8)'],
    products = ['C[CH]CCCC(49)'],
    transitionState = 'TS10',
    kinetics = Arrhenius(A=(4e+13,'cm^3/(mol*s)','*|/',3.16), n=0, Ea=(0,'kJ/mol'), T0=(1,'K'), Tmin=(300,'K'), Tmax=(2000,'K'), comment="""From training reaction 59 used for H_rad;C_rad/H/NonDeC
Exact match found for rate rule [C_rad/H/NonDeC;H_rad]
Euclidian distance = 0
Multiplied by reaction path degeneracy 2.0
family: R_Recombination"""),
)

reaction(
    label = 'reaction26',
    reactants = ['C[CH][CH]CCC(67)', 'H(8)'],
    products = ['C[CH]CCCC(49)'],
    transitionState = 'TS11',
    kinetics = Arrhenius(A=(4.34078e+06,'m^3/(mol*s)'), n=0.278577, Ea=(0,'kJ/mol'), T0=(1,'K'), comment="""Estimated using an average for rate rule [Y_rad;H_rad]
Euclidian distance = 0
family: R_Recombination
Ea raised from -1.4 to 0 kJ/mol."""),
)

reaction(
    label = 'reaction27',
    reactants = ['[CH2]CCC[CH]C(53)', 'H(8)'],
    products = ['C[CH]CCCC(49)'],
    transitionState = 'TS12',
    kinetics = Arrhenius(A=(3.48677e-12,'cm^3/(molecule*s)'), n=0.6, Ea=(0,'kJ/mol'), T0=(1,'K'), comment="""From training reaction 18 used for C_rad/H2/Cs;H_rad
Exact match found for rate rule [C_rad/H2/Cs;H_rad]
Euclidian distance = 0
family: R_Recombination
Ea raised from -3.3 to 0 kJ/mol."""),
)

reaction(
    label = 'reaction28',
    reactants = ['[CH2][CH]CCCC(54)', 'H(8)'],
    products = ['C[CH]CCCC(49)'],
    transitionState = 'TS13',
    kinetics = Arrhenius(A=(4.34078e+06,'m^3/(mol*s)'), n=0.278577, Ea=(0,'kJ/mol'), T0=(1,'K'), comment="""Estimated using an average for rate rule [Y_rad;H_rad]
Euclidian distance = 0
family: R_Recombination
Ea raised from -1.4 to 0 kJ/mol."""),
)

reaction(
    label = 'reaction14',
    reactants = ['CH2(S)(14)', 'C[CH]CCC(25)'],
    products = ['C[CH]CCCC(49)'],
    transitionState = 'TS14',
    kinetics = Arrhenius(A=(1.31021e+06,'m^3/(mol*s)'), n=0.189, Ea=(0,'kJ/mol'), T0=(1,'K'), comment="""Estimated using template [carbene;C_pri] for rate rule [carbene;C_pri/NonDeC]
Euclidian distance = 1.0
Multiplied by reaction path degeneracy 3.0
family: 1,2_Insertion_carbene
Ea raised from -1.5 to 0 kJ/mol."""),
)

reaction(
    label = 'reaction30',
    reactants = ['[CH2]C(C)CCC(68)'],
    products = ['C[CH]CCCC(49)'],
    transitionState = 'TS15',
    kinetics = Arrhenius(A=(6.55606e+10,'s^-1'), n=0.64, Ea=(159.935,'kJ/mol'), T0=(1,'K'), comment="""Estimated using template [cCs(-HC)CJ;CsJ;C] for rate rule [cCs(-HC)CJ;CsJ-HH;C]
Euclidian distance = 1.0
family: 1,2_shiftC"""),
)

reaction(
    label = 'reaction31',
    reactants = ['[CH2]CCC(3)', '[CH]C(20)'],
    products = ['C[CH]CCCC(49)'],
    transitionState = 'TS16',
    kinetics = Arrhenius(A=(759150,'m^3/(mol*s)'), n=0.575199, Ea=(34.3157,'kJ/mol'), T0=(1,'K'), comment="""Estimated using template [Y_rad;Birad] for rate rule [C_rad/H2/Cs;Birad]
Euclidian distance = 3.0
family: Birad_R_Recombination"""),
)

reaction(
    label = 'reaction32',
    reactants = ['[CH3](11)', '[CH]CCCC(33)'],
    products = ['C[CH]CCCC(49)'],
    transitionState = 'TS17',
    kinetics = Arrhenius(A=(759150,'m^3/(mol*s)'), n=0.575199, Ea=(34.3157,'kJ/mol'), T0=(1,'K'), comment="""Estimated using template [Y_rad;Birad] for rate rule [C_methyl;Birad]
Euclidian distance = 2.0
family: Birad_R_Recombination"""),
)

reaction(
    label = 'reaction33',
    reactants = ['C[C]CCCC(69)', 'H(8)'],
    products = ['C[CH]CCCC(49)'],
    transitionState = 'TS18',
    kinetics = Arrhenius(A=(1e+07,'m^3/(mol*s)'), n=0, Ea=(0,'kJ/mol'), T0=(1,'K'), comment="""Estimated using an average for rate rule [H_rad;Birad]
Euclidian distance = 0
family: Birad_R_Recombination"""),
)

reaction(
    label = 'reaction1',
    reactants = ['C=CCCCC(48)', 'H(8)'],
    products = ['[CH2]CCCCC(10)'],
    transitionState = 'TS19',
    kinetics = Arrhenius(A=(3.01e+08,'cm^3/(mol*s)'), n=1.6, Ea=(10.0416,'kJ/mol'), T0=(1,'K'), comment="""From training reaction 10 used for Cds-CsH_Cds-HH;HJ
Exact match found for rate rule [Cds-CsH_Cds-HH;HJ]
Euclidian distance = 0
family: R_Addition_MultipleBond"""),
)

reaction(
    label = 'reaction1',
    reactants = ['[CH2]CCC(3)', 'C=C(17)'],
    products = ['[CH2]CCCCC(10)'],
    transitionState = 'TS20',
    kinetics = Arrhenius(A=(4240,'cm^3/(mol*s)'), n=2.41, Ea=(21.171,'kJ/mol'), T0=(1,'K'), Tmin=(300,'K'), Tmax=(1500,'K'), comment="""From training reaction 220 used for Cds-HH_Cds-HH;CsJ-CsHH
Exact match found for rate rule [Cds-HH_Cds-HH;CsJ-CsHH]
Euclidian distance = 0
Multiplied by reaction path degeneracy 2.0
family: R_Addition_MultipleBond"""),
)

reaction(
    label = 'reaction21',
    reactants = ['[CH2]CCCCC(10)'],
    products = ['CC[CH]CCC(50)'],
    transitionState = 'TS21',
    kinetics = Arrhenius(A=(1.064e+06,'s^-1'), n=1.93, Ea=(141.419,'kJ/mol'), T0=(1,'K'), comment="""From training reaction 108 used for R3H_SS_Cs;C_rad_out_2H;Cs_H_out_H/(NonDeC/Cs)
Exact match found for rate rule [R3H_SS_Cs;C_rad_out_2H;Cs_H_out_H/(NonDeC/Cs)]
Euclidian distance = 0
Multiplied by reaction path degeneracy 2.0
family: intra_H_migration"""),
)

reaction(
    label = 'reaction5',
    reactants = ['[CH2]CC(4)', '[CH2]C[CH2](27)'],
    products = ['[CH2]CCCCC(10)'],
    transitionState = 'TS22',
    kinetics = Arrhenius(A=(1.746e+15,'cm^3/(mol*s)'), n=-0.699, Ea=(0,'kJ/mol'), T0=(1,'K'), Tmin=(200,'K'), Tmax=(2000,'K'), comment="""From training reaction 11 used for C_rad/H2/Cs;C_rad/H2/Cs
Exact match found for rate rule [C_rad/H2/Cs;C_rad/H2/Cs]
Euclidian distance = 0
Multiplied by reaction path degeneracy 2.0
family: R_Recombination
Ea raised from -0.0 to 0 kJ/mol."""),
)

reaction(
    label = 'reaction6',
    reactants = ['[CH2]CC[CH2](29)', 'C[CH2](6)'],
    products = ['[CH2]CCCCC(10)'],
    transitionState = 'TS23',
    kinetics = Arrhenius(A=(1.746e+15,'cm^3/(mol*s)'), n=-0.699, Ea=(0,'kJ/mol'), T0=(1,'K'), Tmin=(200,'K'), Tmax=(2000,'K'), comment="""From training reaction 11 used for C_rad/H2/Cs;C_rad/H2/Cs
Exact match found for rate rule [C_rad/H2/Cs;C_rad/H2/Cs]
Euclidian distance = 0
Multiplied by reaction path degeneracy 2.0
family: R_Recombination
Ea raised from -0.0 to 0 kJ/mol."""),
)

reaction(
    label = 'reaction4',
    reactants = ['[CH2]CC[CH]CC(51)', 'H(8)'],
    products = ['[CH2]CCCCC(10)'],
    transitionState = 'TS24',
    kinetics = Arrhenius(A=(2e+13,'cm^3/(mol*s)','*|/',3.16), n=0, Ea=(0,'kJ/mol'), T0=(1,'K'), Tmin=(300,'K'), Tmax=(2000,'K'), comment="""From training reaction 59 used for H_rad;C_rad/H/NonDeC
Exact match found for rate rule [C_rad/H/NonDeC;H_rad]
Euclidian distance = 0
family: R_Recombination"""),
)

reaction(
    label = 'reaction8',
    reactants = ['[CH2]CCC(3)', '[CH2][CH2](18)'],
    products = ['[CH2]CCCCC(10)'],
    transitionState = 'TS25',
    kinetics = Arrhenius(A=(1.25071e+07,'m^3/(mol*s)'), n=0.093384, Ea=(0.402701,'kJ/mol'), T0=(1,'K'), comment="""Estimated using an average for rate rule [Y_rad;C_rad/H2/Cs]
Euclidian distance = 0
Multiplied by reaction path degeneracy 2.0
family: R_Recombination"""),
)

reaction(
    label = 'reaction5',
    reactants = ['[CH2]C[CH]CCC(52)', 'H(8)'],
    products = ['[CH2]CCCCC(10)'],
    transitionState = 'TS26',
    kinetics = Arrhenius(A=(2e+13,'cm^3/(mol*s)','*|/',3.16), n=0, Ea=(0,'kJ/mol'), T0=(1,'K'), Tmin=(300,'K'), Tmax=(2000,'K'), comment="""From training reaction 59 used for H_rad;C_rad/H/NonDeC
Exact match found for rate rule [C_rad/H/NonDeC;H_rad]
Euclidian distance = 0
family: R_Recombination"""),
)

reaction(
    label = 'reaction6',
    reactants = ['[CH2]CCC[CH2](32)', '[CH3](11)'],
    products = ['[CH2]CCCCC(10)'],
    transitionState = 'TS27',
    kinetics = Arrhenius(A=(2.46e+15,'cm^3/(mol*s)'), n=-0.562, Ea=(0.085772,'kJ/mol'), T0=(1,'K'), Tmin=(200,'K'), Tmax=(2000,'K'), comment="""From training reaction 10 used for C_methyl;C_rad/H2/Cs
Exact match found for rate rule [C_rad/H2/Cs;C_methyl]
Euclidian distance = 0
Multiplied by reaction path degeneracy 2.0
family: R_Recombination"""),
)

reaction(
    label = 'reaction7',
    reactants = ['[CH2]CCC[CH]C(53)', 'H(8)'],
    products = ['[CH2]CCCCC(10)'],
    transitionState = 'TS28',
    kinetics = Arrhenius(A=(2e+13,'cm^3/(mol*s)','*|/',3.16), n=0, Ea=(0,'kJ/mol'), T0=(1,'K'), Tmin=(300,'K'), Tmax=(2000,'K'), comment="""From training reaction 59 used for H_rad;C_rad/H/NonDeC
Exact match found for rate rule [C_rad/H/NonDeC;H_rad]
Euclidian distance = 0
family: R_Recombination"""),
)

reaction(
    label = 'reaction8',
    reactants = ['[CH2][CH]CCCC(54)', 'H(8)'],
    products = ['[CH2]CCCCC(10)'],
    transitionState = 'TS29',
    kinetics = Arrhenius(A=(4.34078e+06,'m^3/(mol*s)'), n=0.278577, Ea=(0,'kJ/mol'), T0=(1,'K'), comment="""Estimated using an average for rate rule [Y_rad;H_rad]
Euclidian distance = 0
family: R_Recombination
Ea raised from -1.4 to 0 kJ/mol."""),
)

reaction(
    label = 'reaction9',
    reactants = ['[CH2]CCCC[CH2](55)', 'H(8)'],
    products = ['[CH2]CCCCC(10)'],
    transitionState = 'TS30',
    kinetics = Arrhenius(A=(6.97354e-12,'cm^3/(molecule*s)'), n=0.6, Ea=(0,'kJ/mol'), T0=(1,'K'), comment="""From training reaction 18 used for C_rad/H2/Cs;H_rad
Exact match found for rate rule [C_rad/H2/Cs;H_rad]
Euclidian distance = 0
Multiplied by reaction path degeneracy 2.0
family: R_Recombination
Ea raised from -3.3 to 0 kJ/mol."""),
)

reaction(
    label = 'reaction31',
    reactants = ['[CH2]CCCC(5)', 'CH2(S)(14)'],
    products = ['[CH2]CCCCC(10)'],
    transitionState = 'TS31',
    kinetics = Arrhenius(A=(1.31021e+06,'m^3/(mol*s)'), n=0.189, Ea=(0,'kJ/mol'), T0=(1,'K'), comment="""Estimated using template [carbene;C_pri] for rate rule [carbene;C_pri/NonDeC]
Euclidian distance = 1.0
Multiplied by reaction path degeneracy 3.0
family: 1,2_Insertion_carbene
Ea raised from -1.5 to 0 kJ/mol."""),
)

reaction(
    label = 'reaction15',
    reactants = ['[CH2]CCCC(5)', 'CH2(T)(19)'],
    products = ['[CH2]CCCCC(10)'],
    transitionState = 'TS32',
    kinetics = Arrhenius(A=(759150,'m^3/(mol*s)'), n=0.575199, Ea=(34.3157,'kJ/mol'), T0=(1,'K'), comment="""Estimated using template [Y_rad;Birad] for rate rule [C_rad/H2/Cs;Birad]
Euclidian distance = 3.0
family: Birad_R_Recombination"""),
)

reaction(
    label = 'reaction10',
    reactants = ['[CH]CCCCC(56)', 'H(8)'],
    products = ['[CH2]CCCCC(10)'],
    transitionState = 'TS33',
    kinetics = Arrhenius(A=(1e+07,'m^3/(mol*s)'), n=0, Ea=(0,'kJ/mol'), T0=(1,'K'), comment="""Estimated using an average for rate rule [H_rad;Birad]
Euclidian distance = 0
family: Birad_R_Recombination"""),
)

reaction(
    label = 'reaction34',
    reactants = ['H(8)', 'CCC=CCC(81)'],
    products = ['CC[CH]CCC(50)'],
    transitionState = 'TS34',
    kinetics = Arrhenius(A=(2.92e+08,'cm^3/(mol*s)'), n=1.64, Ea=(5.73208,'kJ/mol'), T0=(1,'K'), Tmin=(300,'K'), Tmax=(1500,'K'), comment="""From training reaction 2555 used for Cds-CsH_Cds-CsH;HJ
Exact match found for rate rule [Cds-CsH_Cds-CsH;HJ]
Euclidian distance = 0
Multiplied by reaction path degeneracy 2.0
family: R_Addition_MultipleBond"""),
)

reaction(
    label = 'reaction35',
    reactants = ['H(8)', 'CC=CCCC(64)'],
    products = ['CC[CH]CCC(50)'],
    transitionState = 'TS35',
    kinetics = Arrhenius(A=(1.46e+08,'cm^3/(mol*s)'), n=1.64, Ea=(5.73208,'kJ/mol'), T0=(1,'K'), Tmin=(300,'K'), Tmax=(1500,'K'), comment="""From training reaction 2555 used for Cds-CsH_Cds-CsH;HJ
Exact match found for rate rule [Cds-CsH_Cds-CsH;HJ]
Euclidian distance = 0
family: R_Addition_MultipleBond"""),
)

reaction(
    label = 'reaction36',
    reactants = ['C[CH2](6)', 'C=CCC(42)'],
    products = ['CC[CH]CCC(50)'],
    transitionState = 'TS36',
    kinetics = Arrhenius(A=(2130,'cm^3/(mol*s)'), n=2.41, Ea=(19.874,'kJ/mol'), T0=(1,'K'), Tmin=(300,'K'), Tmax=(1500,'K'), comment="""Estimated using an average for rate rule [Cds-HH_Cds-CsH;CsJ-CsHH]
Euclidian distance = 0
family: R_Addition_MultipleBond"""),
)

reaction(
    label = 'reaction37',
    reactants = ['[CH3](11)', 'C=CCCC(24)'],
    products = ['CC[CH]CCC(50)'],
    transitionState = 'TS37',
    kinetics = Arrhenius(A=(0.105698,'m^3/(mol*s)'), n=2.13, Ea=(23.0748,'kJ/mol'), T0=(1,'K'), comment="""Estimated using an average for rate rule [Cds-HH_Cds-CsH;CsJ-HHH]
Euclidian distance = 0
family: R_Addition_MultipleBond"""),
)

reaction(
    label = 'reaction38',
    reactants = ['C[CH2](6)', '[CH2][CH]CC(45)'],
    products = ['CC[CH]CCC(50)'],
    transitionState = 'TS38',
    kinetics = Arrhenius(A=(6.25357e+06,'m^3/(mol*s)'), n=0.093384, Ea=(0.402701,'kJ/mol'), T0=(1,'K'), comment="""Estimated using an average for rate rule [Y_rad;C_rad/H2/Cs]
Euclidian distance = 0
family: R_Recombination"""),
)

reaction(
    label = 'reaction39',
    reactants = ['[CH2]C[CH]CC(28)', '[CH3](11)'],
    products = ['CC[CH]CCC(50)'],
    transitionState = 'TS39',
    kinetics = Arrhenius(A=(1.23e+15,'cm^3/(mol*s)'), n=-0.562, Ea=(0.085772,'kJ/mol'), T0=(1,'K'), Tmin=(200,'K'), Tmax=(2000,'K'), comment="""From training reaction 10 used for C_methyl;C_rad/H2/Cs
Exact match found for rate rule [C_rad/H2/Cs;C_methyl]
Euclidian distance = 0
family: R_Recombination"""),
)

reaction(
    label = 'reaction40',
    reactants = ['H(8)', 'C[CH]C[CH]CC(65)'],
    products = ['CC[CH]CCC(50)'],
    transitionState = 'TS40',
    kinetics = Arrhenius(A=(2e+13,'cm^3/(mol*s)','*|/',3.16), n=0, Ea=(0,'kJ/mol'), T0=(1,'K'), Tmin=(300,'K'), Tmax=(2000,'K'), comment="""From training reaction 59 used for H_rad;C_rad/H/NonDeC
Exact match found for rate rule [C_rad/H/NonDeC;H_rad]
Euclidian distance = 0
family: R_Recombination"""),
)

reaction(
    label = 'reaction41',
    reactants = ['CC[CH][CH]CC(82)', 'H(8)'],
    products = ['CC[CH]CCC(50)'],
    transitionState = 'TS41',
    kinetics = Arrhenius(A=(8.68156e+06,'m^3/(mol*s)'), n=0.278577, Ea=(0,'kJ/mol'), T0=(1,'K'), comment="""Estimated using an average for rate rule [Y_rad;H_rad]
Euclidian distance = 0
Multiplied by reaction path degeneracy 2.0
family: R_Recombination
Ea raised from -1.4 to 0 kJ/mol."""),
)

reaction(
    label = 'reaction42',
    reactants = ['[CH2][CH]CCC(31)', '[CH3](11)'],
    products = ['CC[CH]CCC(50)'],
    transitionState = 'TS42',
    kinetics = Arrhenius(A=(1.66881e+08,'m^3/(mol*s)'), n=-0.401267, Ea=(0,'kJ/mol'), T0=(1,'K'), comment="""Estimated using an average for rate rule [Y_rad;C_methyl]
Euclidian distance = 0
family: R_Recombination
Ea raised from -6.7 to 0 kJ/mol."""),
)

reaction(
    label = 'reaction43',
    reactants = ['C[CH][CH]CCC(67)', 'H(8)'],
    products = ['CC[CH]CCC(50)'],
    transitionState = 'TS43',
    kinetics = Arrhenius(A=(4.34078e+06,'m^3/(mol*s)'), n=0.278577, Ea=(0,'kJ/mol'), T0=(1,'K'), comment="""Estimated using an average for rate rule [Y_rad;H_rad]
Euclidian distance = 0
family: R_Recombination
Ea raised from -1.4 to 0 kJ/mol."""),
)

reaction(
    label = 'reaction44',
    reactants = ['[CH2]CC[CH]CC(51)', 'H(8)'],
    products = ['CC[CH]CCC(50)'],
    transitionState = 'TS44',
    kinetics = Arrhenius(A=(3.48677e-12,'cm^3/(molecule*s)'), n=0.6, Ea=(0,'kJ/mol'), T0=(1,'K'), comment="""From training reaction 18 used for C_rad/H2/Cs;H_rad
Exact match found for rate rule [C_rad/H2/Cs;H_rad]
Euclidian distance = 0
family: R_Recombination
Ea raised from -3.3 to 0 kJ/mol."""),
)

reaction(
    label = 'reaction45',
    reactants = ['[CH2]C[CH]CCC(52)', 'H(8)'],
    products = ['CC[CH]CCC(50)'],
    transitionState = 'TS45',
    kinetics = Arrhenius(A=(3.48677e-12,'cm^3/(molecule*s)'), n=0.6, Ea=(0,'kJ/mol'), T0=(1,'K'), comment="""From training reaction 18 used for C_rad/H2/Cs;H_rad
Exact match found for rate rule [C_rad/H2/Cs;H_rad]
Euclidian distance = 0
family: R_Recombination
Ea raised from -3.3 to 0 kJ/mol."""),
)

reaction(
    label = 'reaction46',
    reactants = ['CH2(S)(14)', 'CC[CH]CC(26)'],
    products = ['CC[CH]CCC(50)'],
    transitionState = 'TS46',
    kinetics = Arrhenius(A=(2.62043e+06,'m^3/(mol*s)'), n=0.189, Ea=(0,'kJ/mol'), T0=(1,'K'), comment="""Estimated using template [carbene;C_pri] for rate rule [carbene;C_pri/NonDeC]
Euclidian distance = 1.0
Multiplied by reaction path degeneracy 6.0
family: 1,2_Insertion_carbene
Ea raised from -1.5 to 0 kJ/mol."""),
)

reaction(
    label = 'reaction47',
    reactants = ['CH2(S)(14)', 'C[CH]CCC(25)'],
    products = ['CC[CH]CCC(50)'],
    transitionState = 'TS47',
    kinetics = Arrhenius(A=(215646,'m^3/(mol*s)'), n=0.444, Ea=(0,'kJ/mol'), T0=(1,'K'), comment="""Estimated using an average for rate rule [carbene;R_H]
Euclidian distance = 0
Multiplied by reaction path degeneracy 3.0
family: 1,2_Insertion_carbene
Ea raised from -5.1 to 0 kJ/mol."""),
)

reaction(
    label = 'reaction48',
    reactants = ['[CH2]C(CC)CC(83)'],
    products = ['CC[CH]CCC(50)'],
    transitionState = 'TS48',
    kinetics = Arrhenius(A=(1.31121e+11,'s^-1'), n=0.64, Ea=(159.935,'kJ/mol'), T0=(1,'K'), comment="""Estimated using template [cCs(-HC)CJ;CsJ;C] for rate rule [cCs(-HC)CJ;CsJ-HH;C]
Euclidian distance = 1.0
Multiplied by reaction path degeneracy 2.0
family: 1,2_shiftC"""),
)

reaction(
    label = 'reaction49',
    reactants = ['[CH2]C(C)CCC(68)'],
    products = ['CC[CH]CCC(50)'],
    transitionState = 'TS49',
    kinetics = Arrhenius(A=(5.59192e+09,'s^-1'), n=1.025, Ea=(194.765,'kJ/mol'), T0=(1,'K'), comment="""Estimated using template [cCs(-HC)CJ;CsJ;CH3] for rate rule [cCs(-HC)CJ;CsJ-HH;CH3]
Euclidian distance = 1.0
family: 1,2_shiftC"""),
)

reaction(
    label = 'reaction50',
    reactants = ['[CH]CC(39)', '[CH2]CC(4)'],
    products = ['CC[CH]CCC(50)'],
    transitionState = 'TS50',
    kinetics = Arrhenius(A=(759150,'m^3/(mol*s)'), n=0.575199, Ea=(34.3157,'kJ/mol'), T0=(1,'K'), comment="""Estimated using template [Y_rad;Birad] for rate rule [C_rad/H2/Cs;Birad]
Euclidian distance = 3.0
family: Birad_R_Recombination"""),
)

reaction(
    label = 'reaction51',
    reactants = ['C[CH2](6)', '[CH]CCC(46)'],
    products = ['CC[CH]CCC(50)'],
    transitionState = 'TS51',
    kinetics = Arrhenius(A=(759150,'m^3/(mol*s)'), n=0.575199, Ea=(34.3157,'kJ/mol'), T0=(1,'K'), comment="""Estimated using template [Y_rad;Birad] for rate rule [C_rad/H2/Cs;Birad]
Euclidian distance = 3.0
family: Birad_R_Recombination"""),
)

reaction(
    label = 'reaction52',
    reactants = ['H(8)', 'CC[C]CCC(84)'],
    products = ['CC[CH]CCC(50)'],
    transitionState = 'TS52',
    kinetics = Arrhenius(A=(1e+07,'m^3/(mol*s)'), n=0, Ea=(0,'kJ/mol'), T0=(1,'K'), comment="""Estimated using an average for rate rule [H_rad;Birad]
Euclidian distance = 0
family: Birad_R_Recombination"""),
)

reaction(
    label = 'reaction53',
    reactants = ['C=C(C)CCC(91)', 'H(8)'],
    products = ['[CH2]C(C)CCC(68)'],
    transitionState = 'TS53',
    kinetics = Arrhenius(A=(0.00507518,'m^3/(mol*s)'), n=2.82235, Ea=(0,'kJ/mol'), T0=(1,'K'), Tmin=(303.03,'K'), Tmax=(2000,'K'), comment="""From training reaction 102 used for Cds-CsCs_Cds-HH;HJ
Exact match found for rate rule [Cds-CsCs_Cds-HH;HJ]
Euclidian distance = 0
family: R_Addition_MultipleBond
Ea raised from -4.8 to 0 kJ/mol."""),
)

reaction(
    label = 'reaction54',
    reactants = ['[CH2]CC(4)', 'C=CC(36)'],
    products = ['[CH2]C(C)CCC(68)'],
    transitionState = 'TS54',
    kinetics = Arrhenius(A=(1020,'cm^3/(mol*s)'), n=2.41, Ea=(27.3634,'kJ/mol'), T0=(1,'K'), Tmin=(300,'K'), Tmax=(1500,'K'), comment="""From training reaction 418 used for Cds-CsH_Cds-HH;CsJ-CsHH
Exact match found for rate rule [Cds-CsH_Cds-HH;CsJ-CsHH]
Euclidian distance = 0
family: R_Addition_MultipleBond"""),
)

reaction(
    label = 'reaction55',
    reactants = ['[CH3](11)', 'C=CCCC(24)'],
    products = ['[CH2]C(C)CCC(68)'],
    transitionState = 'TS55',
    kinetics = Arrhenius(A=(10000,'cm^3/(mol*s)'), n=2.41, Ea=(29.7482,'kJ/mol'), T0=(1,'K'), Tmin=(300,'K'), Tmax=(1500,'K'), comment="""From training reaction 417 used for Cds-CsH_Cds-HH;CsJ-HHH
Exact match found for rate rule [Cds-CsH_Cds-HH;CsJ-HHH]
Euclidian distance = 0
family: R_Addition_MultipleBond"""),
)

reaction(
    label = 'reaction56',
    reactants = ['[CH2]C(C)CCC(68)'],
    products = ['CCC[C](C)C(92)'],
    transitionState = 'TS56',
    kinetics = Arrhenius(A=(5.265e-07,'s^-1'), n=5.639, Ea=(102.68,'kJ/mol'), T0=(1,'K'), comment="""From training reaction 38 used for R2H_S;C_rad_out_2H;Cs_H_out_Cs2
Exact match found for rate rule [R2H_S;C_rad_out_2H;Cs_H_out_Cs2]
Euclidian distance = 0
family: intra_H_migration"""),
)

reaction(
    label = 'reaction57',
    reactants = ['[CH2]C(C)CCC(68)'],
    products = ['CC[CH]C(C)C(93)'],
    transitionState = 'TS57',
    kinetics = Arrhenius(A=(1.064e+06,'s^-1'), n=1.93, Ea=(141.419,'kJ/mol'), T0=(1,'K'), comment="""From training reaction 108 used for R3H_SS_Cs;C_rad_out_2H;Cs_H_out_H/(NonDeC/Cs)
Exact match found for rate rule [R3H_SS_Cs;C_rad_out_2H;Cs_H_out_H/(NonDeC/Cs)]
Euclidian distance = 0
Multiplied by reaction path degeneracy 2.0
family: intra_H_migration"""),
)

reaction(
    label = 'reaction58',
    reactants = ['[CH2]C(C)CCC(68)'],
    products = ['C[CH]CC(C)C(94)'],
    transitionState = 'TS58',
    kinetics = Arrhenius(A=(6.44e+09,'s^-1'), n=0.13, Ea=(86.6088,'kJ/mol'), T0=(1,'K'), Tmin=(300,'K'), Tmax=(1500,'K'), comment="""From training reaction 131 used for R4H_SSS;C_rad_out_2H;Cs_H_out_H/NonDeC
Exact match found for rate rule [R4H_SSS;C_rad_out_2H;Cs_H_out_H/NonDeC]
Euclidian distance = 0
Multiplied by reaction path degeneracy 2.0
family: intra_H_migration"""),
)

reaction(
    label = 'reaction59',
    reactants = ['[CH2]CCC(C)C(95)'],
    products = ['[CH2]C(C)CCC(68)'],
    transitionState = 'TS59',
    kinetics = Arrhenius(A=(137700,'s^-1'), n=1.68, Ea=(52.7184,'kJ/mol'), T0=(1,'K'), comment="""From training reaction 111 used for R5H_CCC;C_rad_out_2H;Cs_H_out_2H
Exact match found for rate rule [R5H_CCC;C_rad_out_2H;Cs_H_out_2H]
Euclidian distance = 0
Multiplied by reaction path degeneracy 6.0
family: intra_H_migration"""),
)

reaction(
    label = 'reaction60',
    reactants = ['[CH2]CC(4)', '[CH2][CH]C(38)'],
    products = ['[CH2]C(C)CCC(68)'],
    transitionState = 'TS60',
    kinetics = Arrhenius(A=(6.25357e+06,'m^3/(mol*s)'), n=0.093384, Ea=(0.402701,'kJ/mol'), T0=(1,'K'), comment="""Estimated using an average for rate rule [Y_rad;C_rad/H2/Cs]
Euclidian distance = 0
family: R_Recombination"""),
)

reaction(
    label = 'reaction61',
    reactants = ['[CH2][CH]CCC(31)', '[CH3](11)'],
    products = ['[CH2]C(C)CCC(68)'],
    transitionState = 'TS61',
    kinetics = Arrhenius(A=(1.66881e+08,'m^3/(mol*s)'), n=-0.401267, Ea=(0,'kJ/mol'), T0=(1,'K'), comment="""Estimated using an average for rate rule [Y_rad;C_methyl]
Euclidian distance = 0
family: R_Recombination
Ea raised from -6.7 to 0 kJ/mol."""),
)

reaction(
    label = 'reaction62',
    reactants = ['[CH2][C](C)CCC(96)', 'H(8)'],
    products = ['[CH2]C(C)CCC(68)'],
    transitionState = 'TS62',
    kinetics = Arrhenius(A=(4.34078e+06,'m^3/(mol*s)'), n=0.278577, Ea=(0,'kJ/mol'), T0=(1,'K'), comment="""Estimated using an average for rate rule [Y_rad;H_rad]
Euclidian distance = 0
family: R_Recombination
Ea raised from -1.4 to 0 kJ/mol."""),
)

reaction(
    label = 'reaction63',
    reactants = ['[CH2]C([CH2])C(97)', 'C[CH2](6)'],
    products = ['[CH2]C(C)CCC(68)'],
    transitionState = 'TS63',
    kinetics = Arrhenius(A=(1.746e+15,'cm^3/(mol*s)'), n=-0.699, Ea=(0,'kJ/mol'), T0=(1,'K'), Tmin=(200,'K'), Tmax=(2000,'K'), comment="""From training reaction 11 used for C_rad/H2/Cs;C_rad/H2/Cs
Exact match found for rate rule [C_rad/H2/Cs;C_rad/H2/Cs]
Euclidian distance = 0
Multiplied by reaction path degeneracy 2.0
family: R_Recombination
Ea raised from -0.0 to 0 kJ/mol."""),
)

reaction(
    label = 'reaction64',
    reactants = ['[CH2]C(C)[CH]CC(98)', 'H(8)'],
    products = ['[CH2]C(C)CCC(68)'],
    transitionState = 'TS64',
    kinetics = Arrhenius(A=(2e+13,'cm^3/(mol*s)','*|/',3.16), n=0, Ea=(0,'kJ/mol'), T0=(1,'K'), Tmin=(300,'K'), Tmax=(2000,'K'), comment="""From training reaction 59 used for H_rad;C_rad/H/NonDeC
Exact match found for rate rule [C_rad/H/NonDeC;H_rad]
Euclidian distance = 0
family: R_Recombination"""),
)

reaction(
    label = 'reaction65',
    reactants = ['[CH2]CC([CH2])C(99)', '[CH3](11)'],
    products = ['[CH2]C(C)CCC(68)'],
    transitionState = 'TS65',
    kinetics = Arrhenius(A=(1.23e+15,'cm^3/(mol*s)'), n=-0.562, Ea=(0.085772,'kJ/mol'), T0=(1,'K'), Tmin=(200,'K'), Tmax=(2000,'K'), comment="""From training reaction 10 used for C_methyl;C_rad/H2/Cs
Exact match found for rate rule [C_rad/H2/Cs;C_methyl]
Euclidian distance = 0
family: R_Recombination"""),
)

reaction(
    label = 'reaction66',
    reactants = ['[CH2]C(C)C[CH]C(100)', 'H(8)'],
    products = ['[CH2]C(C)CCC(68)'],
    transitionState = 'TS66',
    kinetics = Arrhenius(A=(2e+13,'cm^3/(mol*s)','*|/',3.16), n=0, Ea=(0,'kJ/mol'), T0=(1,'K'), Tmin=(300,'K'), Tmax=(2000,'K'), comment="""From training reaction 59 used for H_rad;C_rad/H/NonDeC
Exact match found for rate rule [C_rad/H/NonDeC;H_rad]
Euclidian distance = 0
family: R_Recombination"""),
)

reaction(
    label = 'reaction67',
    reactants = ['[CH2]C([CH2])CCC(101)', 'H(8)'],
    products = ['[CH2]C(C)CCC(68)'],
    transitionState = 'TS67',
    kinetics = Arrhenius(A=(6.97354e-12,'cm^3/(molecule*s)'), n=0.6, Ea=(0,'kJ/mol'), T0=(1,'K'), comment="""From training reaction 18 used for C_rad/H2/Cs;H_rad
Exact match found for rate rule [C_rad/H2/Cs;H_rad]
Euclidian distance = 0
Multiplied by reaction path degeneracy 2.0
family: R_Recombination
Ea raised from -3.3 to 0 kJ/mol."""),
)

reaction(
    label = 'reaction68',
    reactants = ['[CH2]CCC([CH2])C(102)', 'H(8)'],
    products = ['[CH2]C(C)CCC(68)'],
    transitionState = 'TS68',
    kinetics = Arrhenius(A=(3.48677e-12,'cm^3/(molecule*s)'), n=0.6, Ea=(0,'kJ/mol'), T0=(1,'K'), comment="""From training reaction 18 used for C_rad/H2/Cs;H_rad
Exact match found for rate rule [C_rad/H2/Cs;H_rad]
Euclidian distance = 0
family: R_Recombination
Ea raised from -3.3 to 0 kJ/mol."""),
)

reaction(
    label = 'reaction69',
    reactants = ['[CH2]CCCC(5)', 'CH2(S)(14)'],
    products = ['[CH2]C(C)CCC(68)'],
    transitionState = 'TS69',
    kinetics = Arrhenius(A=(143764,'m^3/(mol*s)'), n=0.444, Ea=(0,'kJ/mol'), T0=(1,'K'), comment="""Estimated using an average for rate rule [carbene;R_H]
Euclidian distance = 0
Multiplied by reaction path degeneracy 2.0
family: 1,2_Insertion_carbene
Ea raised from -5.1 to 0 kJ/mol."""),
)

reaction(
    label = 'reaction70',
    reactants = ['[CH2]C(C)CC(103)', 'CH2(S)(14)'],
    products = ['[CH2]C(C)CCC(68)'],
    transitionState = 'TS70',
    kinetics = Arrhenius(A=(1.31021e+06,'m^3/(mol*s)'), n=0.189, Ea=(0,'kJ/mol'), T0=(1,'K'), comment="""Estimated using template [carbene;C_pri] for rate rule [carbene;C_pri/NonDeC]
Euclidian distance = 1.0
Multiplied by reaction path degeneracy 3.0
family: 1,2_Insertion_carbene
Ea raised from -1.5 to 0 kJ/mol."""),
)

reaction(
    label = 'reaction71',
    reactants = ['CH2(T)(19)', 'C[CH]CCC(25)'],
    products = ['[CH2]C(C)CCC(68)'],
    transitionState = 'TS71',
    kinetics = Arrhenius(A=(759150,'m^3/(mol*s)'), n=0.575199, Ea=(34.3157,'kJ/mol'), T0=(1,'K'), comment="""Estimated using template [Y_rad;Birad] for rate rule [C_rad/H/NonDeC;Birad]
Euclidian distance = 3.0
family: Birad_R_Recombination"""),
)

reaction(
    label = 'reaction72',
    reactants = ['[CH]C(C)CCC(104)', 'H(8)'],
    products = ['[CH2]C(C)CCC(68)'],
    transitionState = 'TS72',
    kinetics = Arrhenius(A=(1e+07,'m^3/(mol*s)'), n=0, Ea=(0,'kJ/mol'), T0=(1,'K'), comment="""Estimated using an average for rate rule [H_rad;Birad]
Euclidian distance = 0
family: Birad_R_Recombination"""),
)

reaction(
    label = 'reaction73',
    reactants = ['H(8)', 'CCC=C(C)C(105)'],
    products = ['CCC[C](C)C(92)'],
    transitionState = 'TS73',
    kinetics = Arrhenius(A=(1.69e+08,'cm^3/(mol*s)'), n=1.64, Ea=(3.5564,'kJ/mol'), T0=(1,'K'), Tmin=(300,'K'), Tmax=(1500,'K'), comment="""From training reaction 2556 used for Cds-CsH_Cds-CsCs;HJ
Exact match found for rate rule [Cds-CsH_Cds-CsCs;HJ]
Euclidian distance = 0
family: R_Addition_MultipleBond"""),
)

reaction(
    label = 'reaction74',
    reactants = ['C=C(C)CCC(91)', 'H(8)'],
    products = ['CCC[C](C)C(92)'],
    transitionState = 'TS74',
    kinetics = Arrhenius(A=(2.88e+08,'cm^3/(mol*s)'), n=1.64, Ea=(0,'kJ/mol'), T0=(1,'K'), Tmin=(300,'K'), Tmax=(1500,'K'), comment="""From training reaction 2543 used for Cds-HH_Cds-CsCs;HJ
Exact match found for rate rule [Cds-HH_Cds-CsCs;HJ]
Euclidian distance = 0
family: R_Addition_MultipleBond
Ea raised from -0.0 to 0 kJ/mol."""),
)

reaction(
    label = 'reaction75',
    reactants = ['C=C(C)C(106)', 'C[CH2](6)'],
    products = ['CCC[C](C)C(92)'],
    transitionState = 'TS75',
    kinetics = Arrhenius(A=(2680,'cm^3/(mol*s)'), n=2.41, Ea=(18.2422,'kJ/mol'), T0=(1,'K'), Tmin=(300,'K'), Tmax=(1500,'K'), comment="""Matched reaction 256 C4H8 + C2H5 <=> C6H13-3 in R_Addition_MultipleBond/training
This reaction matched rate rule [Cds-HH_Cds-CsCs;CsJ-CsHH]
family: R_Addition_MultipleBond"""),
)

reaction(
    label = 'reaction76',
    reactants = ['CC[CH]C(C)C(93)'],
    products = ['CCC[C](C)C(92)'],
    transitionState = 'TS76',
    kinetics = Arrhenius(A=(7.25e+10,'s^-1'), n=0.6, Ea=(154.39,'kJ/mol'), T0=(1,'K'), Tmin=(300,'K'), Tmax=(1500,'K'), comment="""From training reaction 150 used for R2H_S;C_rad_out_H/NonDeC;Cs_H_out_Cs2
Exact match found for rate rule [R2H_S;C_rad_out_H/NonDeC;Cs_H_out_Cs2]
Euclidian distance = 0
family: intra_H_migration"""),
)

reaction(
    label = 'reaction77',
    reactants = ['C[CH]CC(C)C(94)'],
    products = ['CCC[C](C)C(92)'],
    transitionState = 'TS77',
    kinetics = Arrhenius(A=(7.27e+09,'s^-1'), n=0.66, Ea=(144.766,'kJ/mol'), T0=(1,'K'), Tmin=(300,'K'), Tmax=(1500,'K'), comment="""Matched reaction 190 C6H13-4 <=> C6H13-5 in intra_H_migration/training
This reaction matched rate rule [R3H_SS_Cs;C_rad_out_Cs2;Cs_H_out_H/NonDeC]
family: intra_H_migration"""),
)

reaction(
    label = 'reaction78',
    reactants = ['[CH2]CCC(C)C(95)'],
    products = ['CCC[C](C)C(92)'],
    transitionState = 'TS78',
    kinetics = Arrhenius(A=(1.86e+10,'s^-1'), n=0.58, Ea=(109.579,'kJ/mol'), T0=(1,'K'), Tmin=(300,'K'), Tmax=(1500,'K'), comment="""Matched reaction 132 C6H13 <=> C6H13-2 in intra_H_migration/training
This reaction matched rate rule [R4H_SSS;C_rad_out_Cs2;Cs_H_out_2H]
family: intra_H_migration"""),
)

reaction(
    label = 'reaction79',
    reactants = ['[CH2][C](C)C(107)', 'C[CH2](6)'],
    products = ['CCC[C](C)C(92)'],
    transitionState = 'TS79',
    kinetics = Arrhenius(A=(6.25357e+06,'m^3/(mol*s)'), n=0.093384, Ea=(0.402701,'kJ/mol'), T0=(1,'K'), comment="""Estimated using an average for rate rule [Y_rad;C_rad/H2/Cs]
Euclidian distance = 0
family: R_Recombination"""),
)

reaction(
    label = 'reaction80',
    reactants = ['[CH2]C[C](C)C(108)', '[CH3](11)'],
    products = ['CCC[C](C)C(92)'],
    transitionState = 'TS80',
    kinetics = Arrhenius(A=(1.23e+15,'cm^3/(mol*s)'), n=-0.562, Ea=(0.085772,'kJ/mol'), T0=(1,'K'), Tmin=(200,'K'), Tmax=(2000,'K'), comment="""From training reaction 10 used for C_methyl;C_rad/H2/Cs
Exact match found for rate rule [C_rad/H2/Cs;C_methyl]
Euclidian distance = 0
family: R_Recombination"""),
)

reaction(
    label = 'reaction81',
    reactants = ['H(8)', 'C[CH]C[C](C)C(109)'],
    products = ['CCC[C](C)C(92)'],
    transitionState = 'TS81',
    kinetics = Arrhenius(A=(2e+13,'cm^3/(mol*s)','*|/',3.16), n=0, Ea=(0,'kJ/mol'), T0=(1,'K'), Tmin=(300,'K'), Tmax=(2000,'K'), comment="""From training reaction 59 used for H_rad;C_rad/H/NonDeC
Exact match found for rate rule [C_rad/H/NonDeC;H_rad]
Euclidian distance = 0
family: R_Recombination"""),
)

reaction(
    label = 'reaction82',
    reactants = ['H(8)', 'CC[CH][C](C)C(110)'],
    products = ['CCC[C](C)C(92)'],
    transitionState = 'TS82',
    kinetics = Arrhenius(A=(4.34078e+06,'m^3/(mol*s)'), n=0.278577, Ea=(0,'kJ/mol'), T0=(1,'K'), comment="""Estimated using an average for rate rule [Y_rad;H_rad]
Euclidian distance = 0
family: R_Recombination
Ea raised from -1.4 to 0 kJ/mol."""),
)

reaction(
    label = 'reaction83',
    reactants = ['H(8)', '[CH2]CC[C](C)C(111)'],
    products = ['CCC[C](C)C(92)'],
    transitionState = 'TS83',
    kinetics = Arrhenius(A=(3.48677e-12,'cm^3/(molecule*s)'), n=0.6, Ea=(0,'kJ/mol'), T0=(1,'K'), comment="""From training reaction 18 used for C_rad/H2/Cs;H_rad
Exact match found for rate rule [C_rad/H2/Cs;H_rad]
Euclidian distance = 0
family: R_Recombination
Ea raised from -3.3 to 0 kJ/mol."""),
)

reaction(
    label = 'reaction84',
    reactants = ['[CH2][C](C)CCC(96)', 'H(8)'],
    products = ['CCC[C](C)C(92)'],
    transitionState = 'TS84',
    kinetics = Arrhenius(A=(4.34078e+06,'m^3/(mol*s)'), n=0.278577, Ea=(0,'kJ/mol'), T0=(1,'K'), comment="""Estimated using an average for rate rule [Y_rad;H_rad]
Euclidian distance = 0
family: R_Recombination
Ea raised from -1.4 to 0 kJ/mol."""),
)

reaction(
    label = 'reaction85',
    reactants = ['CH2(S)(14)', 'CC[C](C)C(112)'],
    products = ['CCC[C](C)C(92)'],
    transitionState = 'TS85',
    kinetics = Arrhenius(A=(1.31021e+06,'m^3/(mol*s)'), n=0.189, Ea=(0,'kJ/mol'), T0=(1,'K'), comment="""Estimated using template [carbene;C_pri] for rate rule [carbene;C_pri/NonDeC]
Euclidian distance = 1.0
Multiplied by reaction path degeneracy 3.0
family: 1,2_Insertion_carbene
Ea raised from -1.5 to 0 kJ/mol."""),
)

reaction(
    label = 'reaction86',
    reactants = ['[CH2]C(C)(C)CC(113)'],
    products = ['CCC[C](C)C(92)'],
    transitionState = 'TS86',
    kinetics = Arrhenius(A=(1.33e+08,'s^-1'), n=1.36, Ea=(157.318,'kJ/mol'), T0=(1,'K'), comment="""Estimated using an average for rate rule [cCs(-R!HR!H)CJ;CsJ-HH;C]
Euclidian distance = 0
family: 1,2_shiftC"""),
)

reaction(
    label = 'reaction87',
    reactants = ['[CH2]CC(4)', 'C[C]C(114)'],
    products = ['CCC[C](C)C(92)'],
    transitionState = 'TS87',
    kinetics = Arrhenius(A=(759150,'m^3/(mol*s)'), n=0.575199, Ea=(34.3157,'kJ/mol'), T0=(1,'K'), comment="""Estimated using template [Y_rad;Birad] for rate rule [C_rad/H2/Cs;Birad]
Euclidian distance = 3.0
family: Birad_R_Recombination"""),
)

reaction(
    label = 'reaction88',
    reactants = ['[CH3](11)', 'C[C]CCC(115)'],
    products = ['CCC[C](C)C(92)'],
    transitionState = 'TS88',
    kinetics = Arrhenius(A=(759150,'m^3/(mol*s)'), n=0.575199, Ea=(34.3157,'kJ/mol'), T0=(1,'K'), comment="""Estimated using template [Y_rad;Birad] for rate rule [C_methyl;Birad]
Euclidian distance = 2.0
family: Birad_R_Recombination"""),
)

reaction(
    label = 'reaction89',
    reactants = ['H(8)', 'C=CCC(C)C(116)'],
    products = ['[CH2]CCC(C)C(95)'],
    transitionState = 'TS89',
    kinetics = Arrhenius(A=(3.01e+08,'cm^3/(mol*s)'), n=1.6, Ea=(10.0416,'kJ/mol'), T0=(1,'K'), comment="""From training reaction 10 used for Cds-CsH_Cds-HH;HJ
Exact match found for rate rule [Cds-CsH_Cds-HH;HJ]
Euclidian distance = 0
family: R_Addition_MultipleBond"""),
)

reaction(
    label = 'reaction90',
    reactants = ['C=C(17)', '[CH2]C(C)C(117)'],
    products = ['[CH2]CCC(C)C(95)'],
    transitionState = 'TS90',
    kinetics = Arrhenius(A=(4240,'cm^3/(mol*s)'), n=2.41, Ea=(21.171,'kJ/mol'), T0=(1,'K'), Tmin=(300,'K'), Tmax=(1500,'K'), comment="""From training reaction 220 used for Cds-HH_Cds-HH;CsJ-CsHH
Exact match found for rate rule [Cds-HH_Cds-HH;CsJ-CsHH]
Euclidian distance = 0
Multiplied by reaction path degeneracy 2.0
family: R_Addition_MultipleBond"""),
)

reaction(
    label = 'reaction91',
    reactants = ['[CH2]CCC(C)C(95)'],
    products = ['C[CH]CC(C)C(94)'],
    transitionState = 'TS91',
    kinetics = Arrhenius(A=(6.48e+07,'s^-1'), n=1.57, Ea=(147.695,'kJ/mol'), T0=(1,'K'), comment="""From training reaction 106 used for R2H_S;C_rad_out_2H;Cs_H_out_H/(NonDeC/Cs)
Exact match found for rate rule [R2H_S;C_rad_out_2H;Cs_H_out_H/(NonDeC/Cs)]
Euclidian distance = 0
Multiplied by reaction path degeneracy 2.0
family: intra_H_migration"""),
)

reaction(
    label = 'reaction92',
    reactants = ['CC[CH]C(C)C(93)'],
    products = ['[CH2]CCC(C)C(95)'],
    transitionState = 'TS92',
    kinetics = Arrhenius(A=(3.57e+08,'s^-1'), n=1.32, Ea=(161.502,'kJ/mol'), T0=(1,'K'), Tmin=(300,'K'), Tmax=(1500,'K'), comment="""From training reaction 187 used for R3H_SS_Cs;C_rad_out_H/NonDeC;Cs_H_out_2H
Exact match found for rate rule [R3H_SS_Cs;C_rad_out_H/NonDeC;Cs_H_out_2H]
Euclidian distance = 0
Multiplied by reaction path degeneracy 3.0
family: intra_H_migration"""),
)

reaction(
    label = 'reaction93',
    reactants = ['C[CH]C(37)', '[CH2]C[CH2](27)'],
    products = ['[CH2]CCC(C)C(95)'],
    transitionState = 'TS93',
    kinetics = Arrhenius(A=(2.3e+14,'cm^3/(mol*s)','*|/',2), n=-0.35, Ea=(0,'kJ/mol'), T0=(1,'K'), Tmin=(300,'K'), Tmax=(2500,'K'), comment="""From training reaction 76 used for C_rad/H2/Cs;C_rad/H/NonDeC
Exact match found for rate rule [C_rad/H2/Cs;C_rad/H/NonDeC]
Euclidian distance = 0
Multiplied by reaction path degeneracy 2.0
family: R_Recombination"""),
)

reaction(
    label = 'reaction94',
    reactants = ['[CH2]CC[CH]C(30)', '[CH3](11)'],
    products = ['[CH2]CCC(C)C(95)'],
    transitionState = 'TS94',
    kinetics = Arrhenius(A=(6.64e+14,'cm^3/(mol*s)'), n=-0.57, Ea=(0,'kJ/mol'), T0=(1,'K'), Tmin=(713,'K'), Tmax=(1800,'K'), comment="""From training reaction 67 used for C_methyl;C_rad/H/NonDeC
Exact match found for rate rule [C_rad/H/NonDeC;C_methyl]
Euclidian distance = 0
family: R_Recombination"""),
)

reaction(
    label = 'reaction95',
    reactants = ['H(8)', '[CH2]CC[C](C)C(111)'],
    products = ['[CH2]CCC(C)C(95)'],
    transitionState = 'TS95',
    kinetics = Arrhenius(A=(6.68468e+06,'m^3/(mol*s)'), n=-0.0135, Ea=(0,'kJ/mol'), T0=(1,'K'), comment="""Estimated using template [C_rad/Cs3;Y_rad] for rate rule [C_rad/Cs3;H_rad]
Euclidian distance = 1.0
family: R_Recombination
Ea raised from -0.9 to 0 kJ/mol."""),
)

reaction(
    label = 'reaction96',
    reactants = ['[CH2][CH2](18)', '[CH2]C(C)C(117)'],
    products = ['[CH2]CCC(C)C(95)'],
    transitionState = 'TS96',
    kinetics = Arrhenius(A=(1.25071e+07,'m^3/(mol*s)'), n=0.093384, Ea=(0.402701,'kJ/mol'), T0=(1,'K'), comment="""Estimated using an average for rate rule [Y_rad;C_rad/H2/Cs]
Euclidian distance = 0
Multiplied by reaction path degeneracy 2.0
family: R_Recombination"""),
)

reaction(
    label = 'reaction97',
    reactants = ['[CH2]C[CH]C(C)C(118)', 'H(8)'],
    products = ['[CH2]CCC(C)C(95)'],
    transitionState = 'TS97',
    kinetics = Arrhenius(A=(2e+13,'cm^3/(mol*s)','*|/',3.16), n=0, Ea=(0,'kJ/mol'), T0=(1,'K'), Tmin=(300,'K'), Tmax=(2000,'K'), comment="""From training reaction 59 used for H_rad;C_rad/H/NonDeC
Exact match found for rate rule [C_rad/H/NonDeC;H_rad]
Euclidian distance = 0
family: R_Recombination"""),
)

reaction(
    label = 'reaction98',
    reactants = ['H(8)', '[CH2][CH]CC(C)C(119)'],
    products = ['[CH2]CCC(C)C(95)'],
    transitionState = 'TS98',
    kinetics = Arrhenius(A=(4.34078e+06,'m^3/(mol*s)'), n=0.278577, Ea=(0,'kJ/mol'), T0=(1,'K'), comment="""Estimated using an average for rate rule [Y_rad;H_rad]
Euclidian distance = 0
family: R_Recombination
Ea raised from -1.4 to 0 kJ/mol."""),
)

reaction(
    label = 'reaction99',
    reactants = ['[CH2]CCC([CH2])C(102)', 'H(8)'],
    products = ['[CH2]CCC(C)C(95)'],
    transitionState = 'TS99',
    kinetics = Arrhenius(A=(3.48677e-12,'cm^3/(molecule*s)'), n=0.6, Ea=(0,'kJ/mol'), T0=(1,'K'), comment="""From training reaction 18 used for C_rad/H2/Cs;H_rad
Exact match found for rate rule [C_rad/H2/Cs;H_rad]
Euclidian distance = 0
family: R_Recombination
Ea raised from -3.3 to 0 kJ/mol."""),
)

reaction(
    label = 'reaction100',
    reactants = ['[CH2]CCCC(5)', 'CH2(S)(14)'],
    products = ['[CH2]CCC(C)C(95)'],
    transitionState = 'TS100',
    kinetics = Arrhenius(A=(873476,'m^3/(mol*s)'), n=0.189, Ea=(0,'kJ/mol'), T0=(1,'K'), comment="""Estimated using template [carbene;Cs_H] for rate rule [carbene;C/H2/NonDeC]
Euclidian distance = 2.0
Multiplied by reaction path degeneracy 2.0
family: 1,2_Insertion_carbene
Ea raised from -1.5 to 0 kJ/mol."""),
)

reaction(
    label = 'reaction101',
    reactants = ['CH2(T)(19)', '[CH2]CC(C)C(120)'],
    products = ['[CH2]CCC(C)C(95)'],
    transitionState = 'TS101',
    kinetics = Arrhenius(A=(759150,'m^3/(mol*s)'), n=0.575199, Ea=(34.3157,'kJ/mol'), T0=(1,'K'), comment="""Estimated using template [Y_rad;Birad] for rate rule [C_rad/H2/Cs;Birad]
Euclidian distance = 3.0
family: Birad_R_Recombination"""),
)

reaction(
    label = 'reaction102',
    reactants = ['H(8)', '[CH]CCC(C)C(121)'],
    products = ['[CH2]CCC(C)C(95)'],
    transitionState = 'TS102',
    kinetics = Arrhenius(A=(1e+07,'m^3/(mol*s)'), n=0, Ea=(0,'kJ/mol'), T0=(1,'K'), comment="""Estimated using an average for rate rule [H_rad;Birad]
Euclidian distance = 0
family: Birad_R_Recombination"""),
)

reaction(
    label = 'reaction103',
    reactants = ['H(8)', 'CCC=C(C)C(105)'],
    products = ['CC[CH]C(C)C(93)'],
    transitionState = 'TS103',
    kinetics = Arrhenius(A=(7.72e+07,'cm^3/(mol*s)'), n=1.64, Ea=(9.07928,'kJ/mol'), T0=(1,'K'), Tmin=(300,'K'), Tmax=(1500,'K'), comment="""From training reaction 2568 used for Cds-CsCs_Cds-CsH;HJ
Exact match found for rate rule [Cds-CsCs_Cds-CsH;HJ]
Euclidian distance = 0
family: R_Addition_MultipleBond"""),
)

reaction(
    label = 'reaction104',
    reactants = ['H(8)', 'CC=CC(C)C(122)'],
    products = ['CC[CH]C(C)C(93)'],
    transitionState = 'TS104',
    kinetics = Arrhenius(A=(1.46e+08,'cm^3/(mol*s)'), n=1.64, Ea=(5.73208,'kJ/mol'), T0=(1,'K'), Tmin=(300,'K'), Tmax=(1500,'K'), comment="""From training reaction 2555 used for Cds-CsH_Cds-CsH;HJ
Exact match found for rate rule [Cds-CsH_Cds-CsH;HJ]
Euclidian distance = 0
family: R_Addition_MultipleBond"""),
)

reaction(
    label = 'reaction105',
    reactants = ['[CH3](11)', 'CC=CCC(123)'],
    products = ['CC[CH]C(C)C(93)'],
    transitionState = 'TS105',
    kinetics = Arrhenius(A=(10100,'cm^3/(mol*s)'), n=2.41, Ea=(28.4512,'kJ/mol'), T0=(1,'K'), Tmin=(300,'K'), Tmax=(1500,'K'), comment="""From training reaction 435 used for Cds-CsH_Cds-CsH;CsJ-HHH
Exact match found for rate rule [Cds-CsH_Cds-CsH;CsJ-HHH]
Euclidian distance = 0
family: R_Addition_MultipleBond"""),
)

reaction(
    label = 'reaction106',
    reactants = ['[CH3](11)', 'C=CC(C)C(124)'],
    products = ['CC[CH]C(C)C(93)'],
    transitionState = 'TS106',
    kinetics = Arrhenius(A=(0.105698,'m^3/(mol*s)'), n=2.13, Ea=(23.0748,'kJ/mol'), T0=(1,'K'), comment="""Estimated using an average for rate rule [Cds-HH_Cds-CsH;CsJ-HHH]
Euclidian distance = 0
family: R_Addition_MultipleBond"""),
)

reaction(
    label = 'reaction107',
    reactants = ['CC[CH]C(C)C(93)'],
    products = ['C[CH]CC(C)C(94)'],
    transitionState = 'TS107',
    kinetics = Arrhenius(A=(6.76e+09,'s^-1'), n=0.88, Ea=(158.992,'kJ/mol'), T0=(1,'K'), Tmin=(300,'K'), Tmax=(1500,'K'), comment="""From training reaction 357 used for R2H_S;C_rad_out_H/NonDeC;Cs_H_out_H/NonDeC
Exact match found for rate rule [R2H_S;C_rad_out_H/NonDeC;Cs_H_out_H/NonDeC]
Euclidian distance = 0
Multiplied by reaction path degeneracy 2.0
family: intra_H_migration"""),
)

reaction(
    label = 'reaction108',
    reactants = ['[CH3](11)', 'C[CH][CH]CC(125)'],
    products = ['CC[CH]C(C)C(93)'],
    transitionState = 'TS108',
    kinetics = Arrhenius(A=(1.66881e+08,'m^3/(mol*s)'), n=-0.401267, Ea=(0,'kJ/mol'), T0=(1,'K'), comment="""Estimated using an average for rate rule [Y_rad;C_methyl]
Euclidian distance = 0
family: R_Recombination
Ea raised from -6.7 to 0 kJ/mol."""),
)

reaction(
    label = 'reaction109',
    reactants = ['H(8)', 'CC[CH][C](C)C(110)'],
    products = ['CC[CH]C(C)C(93)'],
    transitionState = 'TS109',
    kinetics = Arrhenius(A=(4.34078e+06,'m^3/(mol*s)'), n=0.278577, Ea=(0,'kJ/mol'), T0=(1,'K'), comment="""Estimated using an average for rate rule [Y_rad;H_rad]
Euclidian distance = 0
family: R_Recombination
Ea raised from -1.4 to 0 kJ/mol."""),
)

reaction(
    label = 'reaction110',
    reactants = ['[CH3](11)', '[CH2][CH]C(C)C(126)'],
    products = ['CC[CH]C(C)C(93)'],
    transitionState = 'TS110',
    kinetics = Arrhenius(A=(1.66881e+08,'m^3/(mol*s)'), n=-0.401267, Ea=(0,'kJ/mol'), T0=(1,'K'), comment="""Estimated using an average for rate rule [Y_rad;C_methyl]
Euclidian distance = 0
family: R_Recombination
Ea raised from -6.7 to 0 kJ/mol."""),
)

reaction(
    label = 'reaction111',
    reactants = ['H(8)', 'C[CH][CH]C(C)C(127)'],
    products = ['CC[CH]C(C)C(93)'],
    transitionState = 'TS111',
    kinetics = Arrhenius(A=(4.34078e+06,'m^3/(mol*s)'), n=0.278577, Ea=(0,'kJ/mol'), T0=(1,'K'), comment="""Estimated using an average for rate rule [Y_rad;H_rad]
Euclidian distance = 0
family: R_Recombination
Ea raised from -1.4 to 0 kJ/mol."""),
)

reaction(
    label = 'reaction112',
    reactants = ['[CH2]C(C)[CH]CC(98)', 'H(8)'],
    products = ['CC[CH]C(C)C(93)'],
    transitionState = 'TS112',
    kinetics = Arrhenius(A=(3.48677e-12,'cm^3/(molecule*s)'), n=0.6, Ea=(0,'kJ/mol'), T0=(1,'K'), comment="""From training reaction 18 used for C_rad/H2/Cs;H_rad
Exact match found for rate rule [C_rad/H2/Cs;H_rad]
Euclidian distance = 0
family: R_Recombination
Ea raised from -3.3 to 0 kJ/mol."""),
)

reaction(
    label = 'reaction113',
    reactants = ['[CH2]C[CH]C(C)C(118)', 'H(8)'],
    products = ['CC[CH]C(C)C(93)'],
    transitionState = 'TS113',
    kinetics = Arrhenius(A=(3.48677e-12,'cm^3/(molecule*s)'), n=0.6, Ea=(0,'kJ/mol'), T0=(1,'K'), comment="""From training reaction 18 used for C_rad/H2/Cs;H_rad
Exact match found for rate rule [C_rad/H2/Cs;H_rad]
Euclidian distance = 0
family: R_Recombination
Ea raised from -3.3 to 0 kJ/mol."""),
)

reaction(
    label = 'reaction114',
    reactants = ['CH2(S)(14)', 'CC[CH]CC(26)'],
    products = ['CC[CH]C(C)C(93)'],
    transitionState = 'TS114',
    kinetics = Arrhenius(A=(287528,'m^3/(mol*s)'), n=0.444, Ea=(0,'kJ/mol'), T0=(1,'K'), comment="""Estimated using an average for rate rule [carbene;R_H]
Euclidian distance = 0
Multiplied by reaction path degeneracy 4.0
family: 1,2_Insertion_carbene
Ea raised from -5.1 to 0 kJ/mol."""),
)

reaction(
    label = 'reaction115',
    reactants = ['CH2(S)(14)', 'C[CH]C(C)C(128)'],
    products = ['CC[CH]C(C)C(93)'],
    transitionState = 'TS115',
    kinetics = Arrhenius(A=(215646,'m^3/(mol*s)'), n=0.444, Ea=(0,'kJ/mol'), T0=(1,'K'), comment="""Estimated using an average for rate rule [carbene;R_H]
Euclidian distance = 0
Multiplied by reaction path degeneracy 3.0
family: 1,2_Insertion_carbene
Ea raised from -5.1 to 0 kJ/mol."""),
)

reaction(
    label = 'reaction116',
    reactants = ['C[CH]C(C)CC(129)'],
    products = ['CC[CH]C(C)C(93)'],
    transitionState = 'TS116',
    kinetics = Arrhenius(A=(5.59192e+09,'s^-1'), n=1.025, Ea=(194.765,'kJ/mol'), T0=(1,'K'), comment="""Estimated using template [cCs(-HC)CJ;CsJ;CH3] for rate rule [cCs(-HC)CJ;CsJ-CsH;CH3]
Euclidian distance = 1.0
family: 1,2_shiftC"""),
)

reaction(
    label = 'reaction117',
    reactants = ['[CH2]C(C)C(C)C(130)'],
    products = ['CC[CH]C(C)C(93)'],
    transitionState = 'TS117',
    kinetics = Arrhenius(A=(5.59192e+09,'s^-1'), n=1.025, Ea=(194.765,'kJ/mol'), T0=(1,'K'), comment="""Estimated using template [cCs(-HC)CJ;CsJ;CH3] for rate rule [cCs(-HC)CJ;CsJ-HH;CH3]
Euclidian distance = 1.0
family: 1,2_shiftC"""),
)

reaction(
    label = 'reaction118',
    reactants = ['[CH]CC(39)', 'C[CH]C(37)'],
    products = ['CC[CH]C(C)C(93)'],
    transitionState = 'TS118',
    kinetics = Arrhenius(A=(759150,'m^3/(mol*s)'), n=0.575199, Ea=(34.3157,'kJ/mol'), T0=(1,'K'), comment="""Estimated using template [Y_rad;Birad] for rate rule [C_rad/H/NonDeC;Birad]
Euclidian distance = 3.0
family: Birad_R_Recombination"""),
)

reaction(
    label = 'reaction119',
    reactants = ['C[CH2](6)', '[CH]C(C)C(131)'],
    products = ['CC[CH]C(C)C(93)'],
    transitionState = 'TS119',
    kinetics = Arrhenius(A=(759150,'m^3/(mol*s)'), n=0.575199, Ea=(34.3157,'kJ/mol'), T0=(1,'K'), comment="""Estimated using template [Y_rad;Birad] for rate rule [C_rad/H2/Cs;Birad]
Euclidian distance = 3.0
family: Birad_R_Recombination"""),
)

reaction(
    label = 'reaction120',
    reactants = ['H(8)', 'CC[C]C(C)C(132)'],
    products = ['CC[CH]C(C)C(93)'],
    transitionState = 'TS120',
    kinetics = Arrhenius(A=(1e+07,'m^3/(mol*s)'), n=0, Ea=(0,'kJ/mol'), T0=(1,'K'), comment="""Estimated using an average for rate rule [H_rad;Birad]
Euclidian distance = 0
family: Birad_R_Recombination"""),
)

reaction(
    label = 'reaction121',
    reactants = ['C=C(CC)CC(133)', 'H(8)'],
    products = ['[CH2]C(CC)CC(83)'],
    transitionState = 'TS121',
    kinetics = Arrhenius(A=(0.00507518,'m^3/(mol*s)'), n=2.82235, Ea=(0,'kJ/mol'), T0=(1,'K'), Tmin=(303.03,'K'), Tmax=(2000,'K'), comment="""From training reaction 102 used for Cds-CsCs_Cds-HH;HJ
Exact match found for rate rule [Cds-CsCs_Cds-HH;HJ]
Euclidian distance = 0
family: R_Addition_MultipleBond
Ea raised from -4.8 to 0 kJ/mol."""),
)

reaction(
    label = 'reaction122',
    reactants = ['C[CH2](6)', 'C=CCC(42)'],
    products = ['[CH2]C(CC)CC(83)'],
    transitionState = 'TS122',
    kinetics = Arrhenius(A=(1020,'cm^3/(mol*s)'), n=2.41, Ea=(27.3634,'kJ/mol'), T0=(1,'K'), Tmin=(300,'K'), Tmax=(1500,'K'), comment="""From training reaction 418 used for Cds-CsH_Cds-HH;CsJ-CsHH
Exact match found for rate rule [Cds-CsH_Cds-HH;CsJ-CsHH]
Euclidian distance = 0
family: R_Addition_MultipleBond"""),
)

reaction(
    label = 'reaction123',
    reactants = ['[CH2]C(CC)CC(83)'],
    products = ['CC[C](C)CC(134)'],
    transitionState = 'TS123',
    kinetics = Arrhenius(A=(5.265e-07,'s^-1'), n=5.639, Ea=(102.68,'kJ/mol'), T0=(1,'K'), comment="""From training reaction 38 used for R2H_S;C_rad_out_2H;Cs_H_out_Cs2
Exact match found for rate rule [R2H_S;C_rad_out_2H;Cs_H_out_Cs2]
Euclidian distance = 0
family: intra_H_migration"""),
)

reaction(
    label = 'reaction124',
    reactants = ['[CH2]C(CC)CC(83)'],
    products = ['C[CH]C(C)CC(129)'],
    transitionState = 'TS124',
    kinetics = Arrhenius(A=(2.36e+10,'s^-1'), n=0.82, Ea=(146.858,'kJ/mol'), T0=(1,'K'), Tmin=(300,'K'), Tmax=(1500,'K'), comment="""From training reaction 186 used for R3H_SS_Cs;C_rad_out_2H;Cs_H_out_H/NonDeC
Exact match found for rate rule [R3H_SS_Cs;C_rad_out_2H;Cs_H_out_H/NonDeC]
Euclidian distance = 0
Multiplied by reaction path degeneracy 4.0
family: intra_H_migration"""),
)

reaction(
    label = 'reaction125',
    reactants = ['[CH2]C(CC)CC(83)'],
    products = ['[CH2]CC(C)CC(135)'],
    transitionState = 'TS125',
    kinetics = Arrhenius(A=(228000,'s^-1'), n=1.74, Ea=(82.8432,'kJ/mol'), T0=(1,'K'), comment="""From training reaction 109 used for R4H_SSS;C_rad_out_2H;Cs_H_out_2H
Exact match found for rate rule [R4H_SSS;C_rad_out_2H;Cs_H_out_2H]
Euclidian distance = 0
Multiplied by reaction path degeneracy 6.0
family: intra_H_migration"""),
)

reaction(
    label = 'reaction126',
    reactants = ['C[CH2](6)', '[CH2][CH]CC(45)'],
    products = ['[CH2]C(CC)CC(83)'],
    transitionState = 'TS126',
    kinetics = Arrhenius(A=(6.25357e+06,'m^3/(mol*s)'), n=0.093384, Ea=(0.402701,'kJ/mol'), T0=(1,'K'), comment="""Estimated using an average for rate rule [Y_rad;C_rad/H2/Cs]
Euclidian distance = 0
family: R_Recombination"""),
)

reaction(
    label = 'reaction127',
    reactants = ['H(8)', '[CH2][C](CC)CC(136)'],
    products = ['[CH2]C(CC)CC(83)'],
    transitionState = 'TS127',
    kinetics = Arrhenius(A=(4.34078e+06,'m^3/(mol*s)'), n=0.278577, Ea=(0,'kJ/mol'), T0=(1,'K'), comment="""Estimated using an average for rate rule [Y_rad;H_rad]
Euclidian distance = 0
family: R_Recombination
Ea raised from -1.4 to 0 kJ/mol."""),
)

reaction(
    label = 'reaction128',
    reactants = ['[CH3](11)', '[CH2]C([CH2])CC(137)'],
    products = ['[CH2]C(CC)CC(83)'],
    transitionState = 'TS128',
    kinetics = Arrhenius(A=(2.46e+15,'cm^3/(mol*s)'), n=-0.562, Ea=(0.085772,'kJ/mol'), T0=(1,'K'), Tmin=(200,'K'), Tmax=(2000,'K'), comment="""From training reaction 10 used for C_methyl;C_rad/H2/Cs
Exact match found for rate rule [C_rad/H2/Cs;C_methyl]
Euclidian distance = 0
Multiplied by reaction path degeneracy 2.0
family: R_Recombination"""),
)

reaction(
    label = 'reaction129',
    reactants = ['H(8)', '[CH2]C([CH]C)CC(138)'],
    products = ['[CH2]C(CC)CC(83)'],
    transitionState = 'TS129',
    kinetics = Arrhenius(A=(2e+13,'cm^3/(mol*s)','*|/',3.16), n=0, Ea=(0,'kJ/mol'), T0=(1,'K'), Tmin=(300,'K'), Tmax=(2000,'K'), comment="""From training reaction 59 used for H_rad;C_rad/H/NonDeC
Exact match found for rate rule [C_rad/H/NonDeC;H_rad]
Euclidian distance = 0
family: R_Recombination"""),
)

reaction(
    label = 'reaction130',
    reactants = ['[CH2]CC([CH2])CC(139)', 'H(8)'],
    products = ['[CH2]C(CC)CC(83)'],
    transitionState = 'TS130',
    kinetics = Arrhenius(A=(3.48677e-12,'cm^3/(molecule*s)'), n=0.6, Ea=(0,'kJ/mol'), T0=(1,'K'), comment="""From training reaction 18 used for C_rad/H2/Cs;H_rad
Exact match found for rate rule [C_rad/H2/Cs;H_rad]
Euclidian distance = 0
family: R_Recombination
Ea raised from -3.3 to 0 kJ/mol."""),
)

reaction(
    label = 'reaction131',
    reactants = ['[CH2]C(C)CC(103)', 'CH2(S)(14)'],
    products = ['[CH2]C(CC)CC(83)'],
    transitionState = 'TS131',
    kinetics = Arrhenius(A=(1.31021e+06,'m^3/(mol*s)'), n=0.189, Ea=(0,'kJ/mol'), T0=(1,'K'), comment="""Estimated using template [carbene;C_pri] for rate rule [carbene;C_pri/NonDeC]
Euclidian distance = 1.0
Multiplied by reaction path degeneracy 3.0
family: 1,2_Insertion_carbene
Ea raised from -1.5 to 0 kJ/mol."""),
)

reaction(
    label = 'reaction132',
    reactants = ['CH2(T)(19)', 'CC[CH]CC(26)'],
    products = ['[CH2]C(CC)CC(83)'],
    transitionState = 'TS132',
    kinetics = Arrhenius(A=(759150,'m^3/(mol*s)'), n=0.575199, Ea=(34.3157,'kJ/mol'), T0=(1,'K'), comment="""Estimated using template [Y_rad;Birad] for rate rule [C_rad/H/NonDeC;Birad]
Euclidian distance = 3.0
family: Birad_R_Recombination"""),
)

reaction(
    label = 'reaction133',
    reactants = ['[CH]C(CC)CC(140)', 'H(8)'],
    products = ['[CH2]C(CC)CC(83)'],
    transitionState = 'TS133',
    kinetics = Arrhenius(A=(1e+07,'m^3/(mol*s)'), n=0, Ea=(0,'kJ/mol'), T0=(1,'K'), comment="""Estimated using an average for rate rule [H_rad;Birad]
Euclidian distance = 0
family: Birad_R_Recombination"""),
)

reaction(
    label = 'reaction134',
    reactants = ['H(8)', 'CC=CC(C)C(122)'],
    products = ['C[CH]CC(C)C(94)'],
    transitionState = 'TS134',
    kinetics = Arrhenius(A=(1.46e+08,'cm^3/(mol*s)'), n=1.64, Ea=(5.73208,'kJ/mol'), T0=(1,'K'), Tmin=(300,'K'), Tmax=(1500,'K'), comment="""From training reaction 2555 used for Cds-CsH_Cds-CsH;HJ
Exact match found for rate rule [Cds-CsH_Cds-CsH;HJ]
Euclidian distance = 0
family: R_Addition_MultipleBond"""),
)

reaction(
    label = 'reaction135',
    reactants = ['H(8)', 'C=CCC(C)C(116)'],
    products = ['C[CH]CC(C)C(94)'],
    transitionState = 'TS135',
    kinetics = Arrhenius(A=(3.36e+08,'cm^3/(mol*s)'), n=1.56, Ea=(2.5104,'kJ/mol'), T0=(1,'K'), comment="""From training reaction 9 used for Cds-HH_Cds-CsH;HJ
Exact match found for rate rule [Cds-HH_Cds-CsH;HJ]
Euclidian distance = 0
family: R_Addition_MultipleBond"""),
)

reaction(
    label = 'reaction136',
    reactants = ['C=CC(36)', 'C[CH]C(37)'],
    products = ['C[CH]CC(C)C(94)'],
    transitionState = 'TS136',
    kinetics = Arrhenius(A=(1710,'cm^3/(mol*s)'), n=2.41, Ea=(14.8532,'kJ/mol'), T0=(1,'K'), Tmin=(300,'K'), Tmax=(1500,'K'), comment="""Matched reaction 239 propene_1 + C3H7-2 <=> C6H13-2 in R_Addition_MultipleBond/training
This reaction matched rate rule [Cds-HH_Cds-Cs\H3/H;CsJ-CsCsH]
family: R_Addition_MultipleBond"""),
)

reaction(
    label = 'reaction137',
    reactants = ['[CH2][CH]C(38)', 'C[CH]C(37)'],
    products = ['C[CH]CC(C)C(94)'],
    transitionState = 'TS137',
    kinetics = Arrhenius(A=(2.79546e+07,'m^3/(mol*s)'), n=-0.206, Ea=(0,'kJ/mol'), T0=(1,'K'), comment="""Estimated using an average for rate rule [Y_rad;C_rad/H/NonDeC]
Euclidian distance = 0
family: R_Recombination
Ea raised from -0.0 to 0 kJ/mol."""),
)

reaction(
    label = 'reaction138',
    reactants = ['[CH3](11)', 'C[CH]C[CH]C(141)'],
    products = ['C[CH]CC(C)C(94)'],
    transitionState = 'TS138',
    kinetics = Arrhenius(A=(1.328e+15,'cm^3/(mol*s)'), n=-0.57, Ea=(0,'kJ/mol'), T0=(1,'K'), Tmin=(713,'K'), Tmax=(1800,'K'), comment="""From training reaction 67 used for C_methyl;C_rad/H/NonDeC
Exact match found for rate rule [C_rad/H/NonDeC;C_methyl]
Euclidian distance = 0
Multiplied by reaction path degeneracy 2.0
family: R_Recombination"""),
)

reaction(
    label = 'reaction139',
    reactants = ['H(8)', 'C[CH]C[C](C)C(109)'],
    products = ['C[CH]CC(C)C(94)'],
    transitionState = 'TS139',
    kinetics = Arrhenius(A=(6.68468e+06,'m^3/(mol*s)'), n=-0.0135, Ea=(0,'kJ/mol'), T0=(1,'K'), comment="""Estimated using template [C_rad/Cs3;Y_rad] for rate rule [C_rad/Cs3;H_rad]
Euclidian distance = 1.0
family: R_Recombination
Ea raised from -0.9 to 0 kJ/mol."""),
)

reaction(
    label = 'reaction140',
    reactants = ['H(8)', 'C[CH][CH]C(C)C(127)'],
    products = ['C[CH]CC(C)C(94)'],
    transitionState = 'TS140',
    kinetics = Arrhenius(A=(4.34078e+06,'m^3/(mol*s)'), n=0.278577, Ea=(0,'kJ/mol'), T0=(1,'K'), comment="""Estimated using an average for rate rule [Y_rad;H_rad]
Euclidian distance = 0
family: R_Recombination
Ea raised from -1.4 to 0 kJ/mol."""),
)

reaction(
    label = 'reaction141',
    reactants = ['[CH2]C(C)C[CH]C(100)', 'H(8)'],
    products = ['C[CH]CC(C)C(94)'],
    transitionState = 'TS141',
    kinetics = Arrhenius(A=(3.48677e-12,'cm^3/(molecule*s)'), n=0.6, Ea=(0,'kJ/mol'), T0=(1,'K'), comment="""From training reaction 18 used for C_rad/H2/Cs;H_rad
Exact match found for rate rule [C_rad/H2/Cs;H_rad]
Euclidian distance = 0
family: R_Recombination
Ea raised from -3.3 to 0 kJ/mol."""),
)

reaction(
    label = 'reaction142',
    reactants = ['H(8)', '[CH2][CH]CC(C)C(119)'],
    products = ['C[CH]CC(C)C(94)'],
    transitionState = 'TS142',
    kinetics = Arrhenius(A=(4.34078e+06,'m^3/(mol*s)'), n=0.278577, Ea=(0,'kJ/mol'), T0=(1,'K'), comment="""Estimated using an average for rate rule [Y_rad;H_rad]
Euclidian distance = 0
family: R_Recombination
Ea raised from -1.4 to 0 kJ/mol."""),
)

reaction(
    label = 'reaction143',
    reactants = ['CH2(S)(14)', 'C[CH]CCC(25)'],
    products = ['C[CH]CC(C)C(94)'],
    transitionState = 'TS143',
    kinetics = Arrhenius(A=(873476,'m^3/(mol*s)'), n=0.189, Ea=(0,'kJ/mol'), T0=(1,'K'), comment="""Estimated using template [carbene;Cs_H] for rate rule [carbene;C/H2/NonDeC]
Euclidian distance = 2.0
Multiplied by reaction path degeneracy 2.0
family: 1,2_Insertion_carbene
Ea raised from -1.5 to 0 kJ/mol."""),
)

reaction(
    label = 'reaction144',
    reactants = ['[CH2]C(C)C(C)C(130)'],
    products = ['C[CH]CC(C)C(94)'],
    transitionState = 'TS144',
    kinetics = Arrhenius(A=(6.55606e+10,'s^-1'), n=0.64, Ea=(159.935,'kJ/mol'), T0=(1,'K'), comment="""Estimated using template [cCs(-HC)CJ;CsJ;C] for rate rule [cCs(-HC)CJ;CsJ-HH;C]
Euclidian distance = 1.0
family: 1,2_shiftC"""),
)

reaction(
    label = 'reaction145',
    reactants = ['[CH]C(20)', '[CH2]C(C)C(117)'],
    products = ['C[CH]CC(C)C(94)'],
    transitionState = 'TS145',
    kinetics = Arrhenius(A=(759150,'m^3/(mol*s)'), n=0.575199, Ea=(34.3157,'kJ/mol'), T0=(1,'K'), comment="""Estimated using template [Y_rad;Birad] for rate rule [C_rad/H2/Cs;Birad]
Euclidian distance = 3.0
family: Birad_R_Recombination"""),
)

reaction(
    label = 'reaction146',
    reactants = ['[CH3](11)', '[CH]CC(C)C(142)'],
    products = ['C[CH]CC(C)C(94)'],
    transitionState = 'TS146',
    kinetics = Arrhenius(A=(759150,'m^3/(mol*s)'), n=0.575199, Ea=(34.3157,'kJ/mol'), T0=(1,'K'), comment="""Estimated using template [Y_rad;Birad] for rate rule [C_methyl;Birad]
Euclidian distance = 2.0
family: Birad_R_Recombination"""),
)

reaction(
    label = 'reaction147',
    reactants = ['H(8)', 'C[C]CC(C)C(143)'],
    products = ['C[CH]CC(C)C(94)'],
    transitionState = 'TS147',
    kinetics = Arrhenius(A=(1e+07,'m^3/(mol*s)'), n=0, Ea=(0,'kJ/mol'), T0=(1,'K'), comment="""Estimated using an average for rate rule [H_rad;Birad]
Euclidian distance = 0
family: Birad_R_Recombination"""),
)

reaction(
    label = 'reaction148',
    reactants = ['H(8)', 'CC=C(C)CC(144)'],
    products = ['CC[C](C)CC(134)'],
    transitionState = 'TS148',
    kinetics = Arrhenius(A=(1.69e+08,'cm^3/(mol*s)'), n=1.64, Ea=(3.5564,'kJ/mol'), T0=(1,'K'), Tmin=(300,'K'), Tmax=(1500,'K'), comment="""From training reaction 2556 used for Cds-CsH_Cds-CsCs;HJ
Exact match found for rate rule [Cds-CsH_Cds-CsCs;HJ]
Euclidian distance = 0
family: R_Addition_MultipleBond"""),
)

reaction(
    label = 'reaction149',
    reactants = ['C=C(CC)CC(133)', 'H(8)'],
    products = ['CC[C](C)CC(134)'],
    transitionState = 'TS149',
    kinetics = Arrhenius(A=(2.88e+08,'cm^3/(mol*s)'), n=1.64, Ea=(0,'kJ/mol'), T0=(1,'K'), Tmin=(300,'K'), Tmax=(1500,'K'), comment="""From training reaction 2543 used for Cds-HH_Cds-CsCs;HJ
Exact match found for rate rule [Cds-HH_Cds-CsCs;HJ]
Euclidian distance = 0
family: R_Addition_MultipleBond
Ea raised from -0.0 to 0 kJ/mol."""),
)

reaction(
    label = 'reaction150',
    reactants = ['[CH3](11)', 'C=C(C)CC(145)'],
    products = ['CC[C](C)CC(134)'],
    transitionState = 'TS150',
    kinetics = Arrhenius(A=(26400,'cm^3/(mol*s)'), n=2.41, Ea=(20.6271,'kJ/mol'), T0=(1,'K'), Tmin=(300,'K'), Tmax=(1500,'K'), comment="""From training reaction 255 used for Cds-HH_Cds-CsCs;CsJ-HHH
Exact match found for rate rule [Cds-HH_Cds-CsCs;CsJ-HHH]
Euclidian distance = 0
family: R_Addition_MultipleBond"""),
)

reaction(
    label = 'reaction151',
    reactants = ['C[CH]C(C)CC(129)'],
    products = ['CC[C](C)CC(134)'],
    transitionState = 'TS151',
    kinetics = Arrhenius(A=(7.25e+10,'s^-1'), n=0.6, Ea=(154.39,'kJ/mol'), T0=(1,'K'), Tmin=(300,'K'), Tmax=(1500,'K'), comment="""From training reaction 150 used for R2H_S;C_rad_out_H/NonDeC;Cs_H_out_Cs2
Exact match found for rate rule [R2H_S;C_rad_out_H/NonDeC;Cs_H_out_Cs2]
Euclidian distance = 0
family: intra_H_migration"""),
)

reaction(
    label = 'reaction152',
    reactants = ['[CH2]CC(C)CC(135)'],
    products = ['CC[C](C)CC(134)'],
    transitionState = 'TS152',
    kinetics = Arrhenius(A=(2.25e+10,'s^-1'), n=0.66, Ea=(137.654,'kJ/mol'), T0=(1,'K'), Tmin=(300,'K'), Tmax=(1500,'K'), comment="""From training reaction 188 used for R3H_SS_Cs;C_rad_out_2H;Cs_H_out_Cs2
Exact match found for rate rule [R3H_SS_Cs;C_rad_out_2H;Cs_H_out_Cs2]
Euclidian distance = 0
family: intra_H_migration"""),
)

reaction(
    label = 'reaction153',
    reactants = ['[CH3](11)', '[CH2][C](C)CC(146)'],
    products = ['CC[C](C)CC(134)'],
    transitionState = 'TS153',
    kinetics = Arrhenius(A=(1.66881e+08,'m^3/(mol*s)'), n=-0.401267, Ea=(0,'kJ/mol'), T0=(1,'K'), comment="""Estimated using an average for rate rule [Y_rad;C_methyl]
Euclidian distance = 0
family: R_Recombination
Ea raised from -6.7 to 0 kJ/mol."""),
)

reaction(
    label = 'reaction154',
    reactants = ['H(8)', 'C[CH][C](C)CC(147)'],
    products = ['CC[C](C)CC(134)'],
    transitionState = 'TS154',
    kinetics = Arrhenius(A=(4.34078e+06,'m^3/(mol*s)'), n=0.278577, Ea=(0,'kJ/mol'), T0=(1,'K'), comment="""Estimated using an average for rate rule [Y_rad;H_rad]
Euclidian distance = 0
family: R_Recombination
Ea raised from -1.4 to 0 kJ/mol."""),
)

reaction(
    label = 'reaction155',
    reactants = ['H(8)', '[CH2]C[C](C)CC(148)'],
    products = ['CC[C](C)CC(134)'],
    transitionState = 'TS155',
    kinetics = Arrhenius(A=(3.48677e-12,'cm^3/(molecule*s)'), n=0.6, Ea=(0,'kJ/mol'), T0=(1,'K'), comment="""From training reaction 18 used for C_rad/H2/Cs;H_rad
Exact match found for rate rule [C_rad/H2/Cs;H_rad]
Euclidian distance = 0
family: R_Recombination
Ea raised from -3.3 to 0 kJ/mol."""),
)

reaction(
    label = 'reaction156',
    reactants = ['H(8)', '[CH2][C](CC)CC(136)'],
    products = ['CC[C](C)CC(134)'],
    transitionState = 'TS156',
    kinetics = Arrhenius(A=(4.34078e+06,'m^3/(mol*s)'), n=0.278577, Ea=(0,'kJ/mol'), T0=(1,'K'), comment="""Estimated using an average for rate rule [Y_rad;H_rad]
Euclidian distance = 0
family: R_Recombination
Ea raised from -1.4 to 0 kJ/mol."""),
)

reaction(
    label = 'reaction157',
    reactants = ['CH2(S)(14)', 'CC[C](C)C(112)'],
    products = ['CC[C](C)CC(134)'],
    transitionState = 'TS157',
    kinetics = Arrhenius(A=(431291,'m^3/(mol*s)'), n=0.444, Ea=(0,'kJ/mol'), T0=(1,'K'), comment="""Estimated using an average for rate rule [carbene;R_H]
Euclidian distance = 0
Multiplied by reaction path degeneracy 6.0
family: 1,2_Insertion_carbene
Ea raised from -5.1 to 0 kJ/mol."""),
)

reaction(
    label = 'reaction158',
    reactants = ['[CH2]C(C)(C)CC(113)'],
    products = ['CC[C](C)CC(134)'],
    transitionState = 'TS158',
    kinetics = Arrhenius(A=(2.66e+08,'s^-1'), n=1.36, Ea=(157.318,'kJ/mol'), T0=(1,'K'), comment="""From training reaction 2 used for cCs(-R!HR!H)CJ;CsJ-HH;CH3
Exact match found for rate rule [cCs(-R!HR!H)CJ;CsJ-HH;CH3]
Euclidian distance = 0
Multiplied by reaction path degeneracy 2.0
family: 1,2_shiftC"""),
)

reaction(
    label = 'reaction159',
    reactants = ['C[CH2](6)', 'C[C]CC(149)'],
    products = ['CC[C](C)CC(134)'],
    transitionState = 'TS159',
    kinetics = Arrhenius(A=(759150,'m^3/(mol*s)'), n=0.575199, Ea=(34.3157,'kJ/mol'), T0=(1,'K'), comment="""Estimated using template [Y_rad;Birad] for rate rule [C_rad/H2/Cs;Birad]
Euclidian distance = 3.0
family: Birad_R_Recombination"""),
)

reaction(
    label = 'reaction160',
    reactants = ['[CH3](11)', 'CC[C]CC(150)'],
    products = ['CC[C](C)CC(134)'],
    transitionState = 'TS160',
    kinetics = Arrhenius(A=(759150,'m^3/(mol*s)'), n=0.575199, Ea=(34.3157,'kJ/mol'), T0=(1,'K'), comment="""Estimated using template [Y_rad;Birad] for rate rule [C_methyl;Birad]
Euclidian distance = 2.0
family: Birad_R_Recombination"""),
)

reaction(
    label = 'reaction161',
    reactants = ['H(8)', 'CC=C(C)CC(144)'],
    products = ['C[CH]C(C)CC(129)'],
    transitionState = 'TS161',
    kinetics = Arrhenius(A=(7.72e+07,'cm^3/(mol*s)'), n=1.64, Ea=(9.07928,'kJ/mol'), T0=(1,'K'), Tmin=(300,'K'), Tmax=(1500,'K'), comment="""From training reaction 2568 used for Cds-CsCs_Cds-CsH;HJ
Exact match found for rate rule [Cds-CsCs_Cds-CsH;HJ]
Euclidian distance = 0
family: R_Addition_MultipleBond"""),
)

reaction(
    label = 'reaction162',
    reactants = ['H(8)', 'C=CC(C)CC(151)'],
    products = ['C[CH]C(C)CC(129)'],
    transitionState = 'TS162',
    kinetics = Arrhenius(A=(3.36e+08,'cm^3/(mol*s)'), n=1.56, Ea=(2.5104,'kJ/mol'), T0=(1,'K'), comment="""From training reaction 9 used for Cds-HH_Cds-CsH;HJ
Exact match found for rate rule [Cds-HH_Cds-CsH;HJ]
Euclidian distance = 0
family: R_Addition_MultipleBond"""),
)

reaction(
    label = 'reaction163',
    reactants = ['C[CH2](6)', 'CC=CC(152)'],
    products = ['C[CH]C(C)CC(129)'],
    transitionState = 'TS163',
    kinetics = Arrhenius(A=(2040,'cm^3/(mol*s)'), n=2.41, Ea=(26.0663,'kJ/mol'), T0=(1,'K'), Tmin=(300,'K'), Tmax=(1500,'K'), comment="""Matched reaction 436 C2H5 + butene2 <=> C6H13-5 in R_Addition_MultipleBond/training
This reaction matched rate rule [Cds-CsH_Cds-CsH;CsJ-CsHH]
family: R_Addition_MultipleBond"""),
)

reaction(
    label = 'reaction164',
    reactants = ['[CH3](11)', 'CC=CCC(123)'],
    products = ['C[CH]C(C)CC(129)'],
    transitionState = 'TS164',
    kinetics = Arrhenius(A=(10100,'cm^3/(mol*s)'), n=2.41, Ea=(28.4512,'kJ/mol'), T0=(1,'K'), Tmin=(300,'K'), Tmax=(1500,'K'), comment="""From training reaction 435 used for Cds-CsH_Cds-CsH;CsJ-HHH
Exact match found for rate rule [Cds-CsH_Cds-CsH;CsJ-HHH]
Euclidian distance = 0
family: R_Addition_MultipleBond"""),
)

reaction(
    label = 'reaction165',
    reactants = ['C[CH]C(C)CC(129)'],
    products = ['[CH2]CC(C)CC(135)'],
    transitionState = 'TS165',
    kinetics = Arrhenius(A=(5.04e-11,'s^-1'), n=6.833, Ea=(117.248,'kJ/mol'), T0=(1,'K'), comment="""From training reaction 33 used for R2H_S;C_rad_out_H/NonDeC;Cs_H_out_2H
Exact match found for rate rule [R2H_S;C_rad_out_H/NonDeC;Cs_H_out_2H]
Euclidian distance = 0
Multiplied by reaction path degeneracy 3.0
family: intra_H_migration"""),
)

reaction(
    label = 'reaction166',
    reactants = ['C[CH2](6)', 'C[CH][CH]C(153)'],
    products = ['C[CH]C(C)CC(129)'],
    transitionState = 'TS166',
    kinetics = Arrhenius(A=(1.25071e+07,'m^3/(mol*s)'), n=0.093384, Ea=(0.402701,'kJ/mol'), T0=(1,'K'), comment="""Estimated using an average for rate rule [Y_rad;C_rad/H2/Cs]
Euclidian distance = 0
Multiplied by reaction path degeneracy 2.0
family: R_Recombination"""),
)

reaction(
    label = 'reaction167',
    reactants = ['[CH3](11)', 'C[CH][CH]CC(125)'],
    products = ['C[CH]C(C)CC(129)'],
    transitionState = 'TS167',
    kinetics = Arrhenius(A=(1.66881e+08,'m^3/(mol*s)'), n=-0.401267, Ea=(0,'kJ/mol'), T0=(1,'K'), comment="""Estimated using an average for rate rule [Y_rad;C_methyl]
Euclidian distance = 0
family: R_Recombination
Ea raised from -6.7 to 0 kJ/mol."""),
)

reaction(
    label = 'reaction168',
    reactants = ['H(8)', 'C[CH][C](C)CC(147)'],
    products = ['C[CH]C(C)CC(129)'],
    transitionState = 'TS168',
    kinetics = Arrhenius(A=(4.34078e+06,'m^3/(mol*s)'), n=0.278577, Ea=(0,'kJ/mol'), T0=(1,'K'), comment="""Estimated using an average for rate rule [Y_rad;H_rad]
Euclidian distance = 0
family: R_Recombination
Ea raised from -1.4 to 0 kJ/mol."""),
)

reaction(
    label = 'reaction169',
    reactants = ['[CH3](11)', '[CH2]C(C)[CH]C(154)'],
    products = ['C[CH]C(C)CC(129)'],
    transitionState = 'TS169',
    kinetics = Arrhenius(A=(1.23e+15,'cm^3/(mol*s)'), n=-0.562, Ea=(0.085772,'kJ/mol'), T0=(1,'K'), Tmin=(200,'K'), Tmax=(2000,'K'), comment="""From training reaction 10 used for C_methyl;C_rad/H2/Cs
Exact match found for rate rule [C_rad/H2/Cs;C_methyl]
Euclidian distance = 0
family: R_Recombination"""),
)

reaction(
    label = 'reaction170',
    reactants = ['H(8)', 'C[CH]C(C)[CH]C(155)'],
    products = ['C[CH]C(C)CC(129)'],
    transitionState = 'TS170',
    kinetics = Arrhenius(A=(4e+13,'cm^3/(mol*s)','*|/',3.16), n=0, Ea=(0,'kJ/mol'), T0=(1,'K'), Tmin=(300,'K'), Tmax=(2000,'K'), comment="""From training reaction 59 used for H_rad;C_rad/H/NonDeC
Exact match found for rate rule [C_rad/H/NonDeC;H_rad]
Euclidian distance = 0
Multiplied by reaction path degeneracy 2.0
family: R_Recombination"""),
)

reaction(
    label = 'reaction171',
    reactants = ['H(8)', '[CH2]C([CH]C)CC(138)'],
    products = ['C[CH]C(C)CC(129)'],
    transitionState = 'TS171',
    kinetics = Arrhenius(A=(3.48677e-12,'cm^3/(molecule*s)'), n=0.6, Ea=(0,'kJ/mol'), T0=(1,'K'), comment="""From training reaction 18 used for C_rad/H2/Cs;H_rad
Exact match found for rate rule [C_rad/H2/Cs;H_rad]
Euclidian distance = 0
family: R_Recombination
Ea raised from -3.3 to 0 kJ/mol."""),
)

reaction(
    label = 'reaction172',
    reactants = ['H(8)', '[CH2]CC(C)[CH]C(156)'],
    products = ['C[CH]C(C)CC(129)'],
    transitionState = 'TS172',
    kinetics = Arrhenius(A=(3.48677e-12,'cm^3/(molecule*s)'), n=0.6, Ea=(0,'kJ/mol'), T0=(1,'K'), comment="""From training reaction 18 used for C_rad/H2/Cs;H_rad
Exact match found for rate rule [C_rad/H2/Cs;H_rad]
Euclidian distance = 0
family: R_Recombination
Ea raised from -3.3 to 0 kJ/mol."""),
)

reaction(
    label = 'reaction173',
    reactants = ['H(8)', '[CH2][CH]C(C)CC(157)'],
    products = ['C[CH]C(C)CC(129)'],
    transitionState = 'TS173',
    kinetics = Arrhenius(A=(4.34078e+06,'m^3/(mol*s)'), n=0.278577, Ea=(0,'kJ/mol'), T0=(1,'K'), comment="""Estimated using an average for rate rule [Y_rad;H_rad]
Euclidian distance = 0
family: R_Recombination
Ea raised from -1.4 to 0 kJ/mol."""),
)

reaction(
    label = 'reaction174',
    reactants = ['CH2(S)(14)', 'C[CH]CCC(25)'],
    products = ['C[CH]C(C)CC(129)'],
    transitionState = 'TS174',
    kinetics = Arrhenius(A=(143764,'m^3/(mol*s)'), n=0.444, Ea=(0,'kJ/mol'), T0=(1,'K'), comment="""Estimated using an average for rate rule [carbene;R_H]
Euclidian distance = 0
Multiplied by reaction path degeneracy 2.0
family: 1,2_Insertion_carbene
Ea raised from -5.1 to 0 kJ/mol."""),
)

reaction(
    label = 'reaction175',
    reactants = ['CH2(S)(14)', 'C[CH]C(C)C(128)'],
    products = ['C[CH]C(C)CC(129)'],
    transitionState = 'TS175',
    kinetics = Arrhenius(A=(2.62043e+06,'m^3/(mol*s)'), n=0.189, Ea=(0,'kJ/mol'), T0=(1,'K'), comment="""Estimated using template [carbene;C_pri] for rate rule [carbene;C_pri/NonDeC]
Euclidian distance = 1.0
Multiplied by reaction path degeneracy 6.0
family: 1,2_Insertion_carbene
Ea raised from -1.5 to 0 kJ/mol."""),
)

reaction(
    label = 'reaction176',
    reactants = ['C[CH]CC(43)', '[CH]C(20)'],
    products = ['C[CH]C(C)CC(129)'],
    transitionState = 'TS176',
    kinetics = Arrhenius(A=(759150,'m^3/(mol*s)'), n=0.575199, Ea=(34.3157,'kJ/mol'), T0=(1,'K'), comment="""Estimated using template [Y_rad;Birad] for rate rule [C_rad/H/NonDeC;Birad]
Euclidian distance = 3.0
family: Birad_R_Recombination"""),
)

reaction(
    label = 'reaction177',
    reactants = ['[CH3](11)', '[CH]C(C)CC(158)'],
    products = ['C[CH]C(C)CC(129)'],
    transitionState = 'TS177',
    kinetics = Arrhenius(A=(759150,'m^3/(mol*s)'), n=0.575199, Ea=(34.3157,'kJ/mol'), T0=(1,'K'), comment="""Estimated using template [Y_rad;Birad] for rate rule [C_methyl;Birad]
Euclidian distance = 2.0
family: Birad_R_Recombination"""),
)

reaction(
    label = 'reaction178',
    reactants = ['H(8)', 'C[C]C(C)CC(159)'],
    products = ['C[CH]C(C)CC(129)'],
    transitionState = 'TS178',
    kinetics = Arrhenius(A=(1e+07,'m^3/(mol*s)'), n=0, Ea=(0,'kJ/mol'), T0=(1,'K'), comment="""Estimated using an average for rate rule [H_rad;Birad]
Euclidian distance = 0
family: Birad_R_Recombination"""),
)

reaction(
    label = 'reaction179',
    reactants = ['C=C(C)C(106)', 'C[CH2](6)'],
    products = ['[CH2]C(C)(C)CC(113)'],
    transitionState = 'TS179',
    kinetics = Arrhenius(A=(634,'cm^3/(mol*s)'), n=2.41, Ea=(31.2545,'kJ/mol'), T0=(1,'K'), Tmin=(300,'K'), Tmax=(1500,'K'), comment="""Matched reaction 616 C4H8-2 + C2H5 <=> C6H13-7 in R_Addition_MultipleBond/training
This reaction matched rate rule [Cds-CsCs_Cds-HH;CsJ-CsHH]
family: R_Addition_MultipleBond"""),
)

reaction(
    label = 'reaction180',
    reactants = ['[CH3](11)', 'C=C(C)CC(145)'],
    products = ['[CH2]C(C)(C)CC(113)'],
    transitionState = 'TS180',
    kinetics = Arrhenius(A=(6240,'cm^3/(mol*s)'), n=2.41, Ea=(33.6394,'kJ/mol'), T0=(1,'K'), Tmin=(300,'K'), Tmax=(1500,'K'), comment="""From training reaction 615 used for Cds-CsCs_Cds-HH;CsJ-HHH
Exact match found for rate rule [Cds-CsCs_Cds-HH;CsJ-HHH]
Euclidian distance = 0
family: R_Addition_MultipleBond"""),
)

reaction(
    label = 'reaction181',
    reactants = ['[CH2]C(C)(C)CC(113)'],
    products = ['C[CH]C(C)(C)C(160)'],
    transitionState = 'TS181',
    kinetics = Arrhenius(A=(1.18e+10,'s^-1'), n=0.82, Ea=(146.858,'kJ/mol'), T0=(1,'K'), Tmin=(300,'K'), Tmax=(1500,'K'), comment="""From training reaction 186 used for R3H_SS_Cs;C_rad_out_2H;Cs_H_out_H/NonDeC
Exact match found for rate rule [R3H_SS_Cs;C_rad_out_2H;Cs_H_out_H/NonDeC]
Euclidian distance = 0
Multiplied by reaction path degeneracy 2.0
family: intra_H_migration"""),
)

reaction(
    label = 'reaction182',
    reactants = ['[CH2]CC(C)(C)C(161)'],
    products = ['[CH2]C(C)(C)CC(113)'],
    transitionState = 'TS182',
    kinetics = Arrhenius(A=(342000,'s^-1'), n=1.74, Ea=(82.8432,'kJ/mol'), T0=(1,'K'), comment="""From training reaction 109 used for R4H_SSS;C_rad_out_2H;Cs_H_out_2H
Exact match found for rate rule [R4H_SSS;C_rad_out_2H;Cs_H_out_2H]
Euclidian distance = 0
Multiplied by reaction path degeneracy 9.0
family: intra_H_migration"""),
)

reaction(
    label = 'reaction183',
    reactants = ['[CH2][C](C)C(107)', 'C[CH2](6)'],
    products = ['[CH2]C(C)(C)CC(113)'],
    transitionState = 'TS183',
    kinetics = Arrhenius(A=(6.25357e+06,'m^3/(mol*s)'), n=0.093384, Ea=(0.402701,'kJ/mol'), T0=(1,'K'), comment="""Estimated using an average for rate rule [Y_rad;C_rad/H2/Cs]
Euclidian distance = 0
family: R_Recombination"""),
)

reaction(
    label = 'reaction184',
    reactants = ['[CH3](11)', '[CH2][C](C)CC(146)'],
    products = ['[CH2]C(C)(C)CC(113)'],
    transitionState = 'TS184',
    kinetics = Arrhenius(A=(1.66881e+08,'m^3/(mol*s)'), n=-0.401267, Ea=(0,'kJ/mol'), T0=(1,'K'), comment="""Estimated using an average for rate rule [Y_rad;C_methyl]
Euclidian distance = 0
family: R_Recombination
Ea raised from -6.7 to 0 kJ/mol."""),
)

reaction(
    label = 'reaction185',
    reactants = ['[CH3](11)', '[CH2]C([CH2])(C)C(162)'],
    products = ['[CH2]C(C)(C)CC(113)'],
    transitionState = 'TS185',
    kinetics = Arrhenius(A=(2.46e+15,'cm^3/(mol*s)'), n=-0.562, Ea=(0.085772,'kJ/mol'), T0=(1,'K'), Tmin=(200,'K'), Tmax=(2000,'K'), comment="""From training reaction 10 used for C_methyl;C_rad/H2/Cs
Exact match found for rate rule [C_rad/H2/Cs;C_methyl]
Euclidian distance = 0
Multiplied by reaction path degeneracy 2.0
family: R_Recombination"""),
)

reaction(
    label = 'reaction186',
    reactants = ['H(8)', '[CH2]C(C)(C)[CH]C(163)'],
    products = ['[CH2]C(C)(C)CC(113)'],
    transitionState = 'TS186',
    kinetics = Arrhenius(A=(2e+13,'cm^3/(mol*s)','*|/',3.16), n=0, Ea=(0,'kJ/mol'), T0=(1,'K'), Tmin=(300,'K'), Tmax=(2000,'K'), comment="""From training reaction 59 used for H_rad;C_rad/H/NonDeC
Exact match found for rate rule [C_rad/H/NonDeC;H_rad]
Euclidian distance = 0
family: R_Recombination"""),
)

reaction(
    label = 'reaction187',
    reactants = ['H(8)', '[CH2]C([CH2])(C)CC(164)'],
    products = ['[CH2]C(C)(C)CC(113)'],
    transitionState = 'TS187',
    kinetics = Arrhenius(A=(6.97354e-12,'cm^3/(molecule*s)'), n=0.6, Ea=(0,'kJ/mol'), T0=(1,'K'), comment="""From training reaction 18 used for C_rad/H2/Cs;H_rad
Exact match found for rate rule [C_rad/H2/Cs;H_rad]
Euclidian distance = 0
Multiplied by reaction path degeneracy 2.0
family: R_Recombination
Ea raised from -3.3 to 0 kJ/mol."""),
)

reaction(
    label = 'reaction188',
    reactants = ['H(8)', '[CH2]CC([CH2])(C)C(165)'],
    products = ['[CH2]C(C)(C)CC(113)'],
    transitionState = 'TS188',
    kinetics = Arrhenius(A=(3.48677e-12,'cm^3/(molecule*s)'), n=0.6, Ea=(0,'kJ/mol'), T0=(1,'K'), comment="""From training reaction 18 used for C_rad/H2/Cs;H_rad
Exact match found for rate rule [C_rad/H2/Cs;H_rad]
Euclidian distance = 0
family: R_Recombination
Ea raised from -3.3 to 0 kJ/mol."""),
)

reaction(
    label = 'reaction189',
    reactants = ['[CH2]C(C)CC(103)', 'CH2(S)(14)'],
    products = ['[CH2]C(C)(C)CC(113)'],
    transitionState = 'TS189',
    kinetics = Arrhenius(A=(71881.9,'m^3/(mol*s)'), n=0.444, Ea=(0,'kJ/mol'), T0=(1,'K'), comment="""Estimated using an average for rate rule [carbene;R_H]
Euclidian distance = 0
family: 1,2_Insertion_carbene
Ea raised from -5.1 to 0 kJ/mol."""),
)

reaction(
    label = 'reaction190',
    reactants = ['CH2(S)(14)', '[CH2]C(C)(C)C(166)'],
    products = ['[CH2]C(C)(C)CC(113)'],
    transitionState = 'TS190',
    kinetics = Arrhenius(A=(3.93064e+06,'m^3/(mol*s)'), n=0.189, Ea=(0,'kJ/mol'), T0=(1,'K'), comment="""Estimated using template [carbene;C_pri] for rate rule [carbene;C_pri/NonDeC]
Euclidian distance = 1.0
Multiplied by reaction path degeneracy 9.0
family: 1,2_Insertion_carbene
Ea raised from -1.5 to 0 kJ/mol."""),
)

reaction(
    label = 'reaction191',
    reactants = ['CH2(T)(19)', 'CC[C](C)C(112)'],
    products = ['[CH2]C(C)(C)CC(113)'],
    transitionState = 'TS191',
    kinetics = Arrhenius(A=(759150,'m^3/(mol*s)'), n=0.575199, Ea=(34.3157,'kJ/mol'), T0=(1,'K'), comment="""Estimated using template [Y_rad;Birad] for rate rule [C_rad/Cs3;Birad]
Euclidian distance = 4.0
family: Birad_R_Recombination"""),
)

reaction(
    label = 'reaction192',
    reactants = ['H(8)', '[CH]C(C)(C)CC(167)'],
    products = ['[CH2]C(C)(C)CC(113)'],
    transitionState = 'TS192',
    kinetics = Arrhenius(A=(1e+07,'m^3/(mol*s)'), n=0, Ea=(0,'kJ/mol'), T0=(1,'K'), comment="""Estimated using an average for rate rule [H_rad;Birad]
Euclidian distance = 0
family: Birad_R_Recombination"""),
)

reaction(
    label = 'reaction193',
    reactants = ['H(8)', 'C=CC(C)CC(151)'],
    products = ['[CH2]CC(C)CC(135)'],
    transitionState = 'TS193',
    kinetics = Arrhenius(A=(3.01e+08,'cm^3/(mol*s)'), n=1.6, Ea=(10.0416,'kJ/mol'), T0=(1,'K'), comment="""From training reaction 10 used for Cds-CsH_Cds-HH;HJ
Exact match found for rate rule [Cds-CsH_Cds-HH;HJ]
Euclidian distance = 0
family: R_Addition_MultipleBond"""),
)

reaction(
    label = 'reaction194',
    reactants = ['C[CH]CC(43)', 'C=C(17)'],
    products = ['[CH2]CC(C)CC(135)'],
    transitionState = 'TS194',
    kinetics = Arrhenius(A=(3400,'cm^3/(mol*s)'), n=2.41, Ea=(16.1921,'kJ/mol'), T0=(1,'K'), Tmin=(300,'K'), Tmax=(1500,'K'), comment="""From training reaction 221 used for Cds-HH_Cds-HH;CsJ-CsCsH
Exact match found for rate rule [Cds-HH_Cds-HH;CsJ-CsCsH]
Euclidian distance = 0
Multiplied by reaction path degeneracy 2.0
family: R_Addition_MultipleBond"""),
)

reaction(
    label = 'reaction195',
    reactants = ['[CH2]C[CH]C(44)', 'C[CH2](6)'],
    products = ['[CH2]CC(C)CC(135)'],
    transitionState = 'TS195',
    kinetics = Arrhenius(A=(1.15e+14,'cm^3/(mol*s)','*|/',2), n=-0.35, Ea=(0,'kJ/mol'), T0=(1,'K'), Tmin=(300,'K'), Tmax=(2500,'K'), comment="""From training reaction 76 used for C_rad/H2/Cs;C_rad/H/NonDeC
Exact match found for rate rule [C_rad/H/NonDeC;C_rad/H2/Cs]
Euclidian distance = 0
family: R_Recombination"""),
)

reaction(
    label = 'reaction196',
    reactants = ['C[CH]CC(43)', '[CH2][CH2](18)'],
    products = ['[CH2]CC(C)CC(135)'],
    transitionState = 'TS196',
    kinetics = Arrhenius(A=(5.59093e+07,'m^3/(mol*s)'), n=-0.206, Ea=(0,'kJ/mol'), T0=(1,'K'), comment="""Estimated using an average for rate rule [Y_rad;C_rad/H/NonDeC]
Euclidian distance = 0
Multiplied by reaction path degeneracy 2.0
family: R_Recombination
Ea raised from -0.0 to 0 kJ/mol."""),
)

reaction(
    label = 'reaction197',
    reactants = ['[CH2]C[CH]CC(28)', '[CH3](11)'],
    products = ['[CH2]CC(C)CC(135)'],
    transitionState = 'TS197',
    kinetics = Arrhenius(A=(6.64e+14,'cm^3/(mol*s)'), n=-0.57, Ea=(0,'kJ/mol'), T0=(1,'K'), Tmin=(713,'K'), Tmax=(1800,'K'), comment="""From training reaction 67 used for C_methyl;C_rad/H/NonDeC
Exact match found for rate rule [C_rad/H/NonDeC;C_methyl]
Euclidian distance = 0
family: R_Recombination"""),
)

reaction(
    label = 'reaction198',
    reactants = ['H(8)', '[CH2]C[C](C)CC(148)'],
    products = ['[CH2]CC(C)CC(135)'],
    transitionState = 'TS198',
    kinetics = Arrhenius(A=(6.68468e+06,'m^3/(mol*s)'), n=-0.0135, Ea=(0,'kJ/mol'), T0=(1,'K'), comment="""Estimated using template [C_rad/Cs3;Y_rad] for rate rule [C_rad/Cs3;H_rad]
Euclidian distance = 1.0
family: R_Recombination
Ea raised from -0.9 to 0 kJ/mol."""),
)

reaction(
    label = 'reaction199',
    reactants = ['[CH2]CC([CH2])C(99)', '[CH3](11)'],
    products = ['[CH2]CC(C)CC(135)'],
    transitionState = 'TS199',
    kinetics = Arrhenius(A=(1.23e+15,'cm^3/(mol*s)'), n=-0.562, Ea=(0.085772,'kJ/mol'), T0=(1,'K'), Tmin=(200,'K'), Tmax=(2000,'K'), comment="""From training reaction 10 used for C_methyl;C_rad/H2/Cs
Exact match found for rate rule [C_rad/H2/Cs;C_methyl]
Euclidian distance = 0
family: R_Recombination"""),
)

reaction(
    label = 'reaction200',
    reactants = ['H(8)', '[CH2]CC(C)[CH]C(156)'],
    products = ['[CH2]CC(C)CC(135)'],
    transitionState = 'TS200',
    kinetics = Arrhenius(A=(2e+13,'cm^3/(mol*s)','*|/',3.16), n=0, Ea=(0,'kJ/mol'), T0=(1,'K'), Tmin=(300,'K'), Tmax=(2000,'K'), comment="""From training reaction 59 used for H_rad;C_rad/H/NonDeC
Exact match found for rate rule [C_rad/H/NonDeC;H_rad]
Euclidian distance = 0
family: R_Recombination"""),
)

reaction(
    label = 'reaction201',
    reactants = ['H(8)', '[CH2][CH]C(C)CC(157)'],
    products = ['[CH2]CC(C)CC(135)'],
    transitionState = 'TS201',
    kinetics = Arrhenius(A=(4.34078e+06,'m^3/(mol*s)'), n=0.278577, Ea=(0,'kJ/mol'), T0=(1,'K'), comment="""Estimated using an average for rate rule [Y_rad;H_rad]
Euclidian distance = 0
family: R_Recombination
Ea raised from -1.4 to 0 kJ/mol."""),
)

reaction(
    label = 'reaction202',
    reactants = ['[CH2]CC([CH2])CC(139)', 'H(8)'],
    products = ['[CH2]CC(C)CC(135)'],
    transitionState = 'TS202',
    kinetics = Arrhenius(A=(3.48677e-12,'cm^3/(molecule*s)'), n=0.6, Ea=(0,'kJ/mol'), T0=(1,'K'), comment="""From training reaction 18 used for C_rad/H2/Cs;H_rad
Exact match found for rate rule [C_rad/H2/Cs;H_rad]
Euclidian distance = 0
family: R_Recombination
Ea raised from -3.3 to 0 kJ/mol."""),
)

reaction(
    label = 'reaction203',
    reactants = ['[CH2]CC(C)C[CH2](168)', 'H(8)'],
    products = ['[CH2]CC(C)CC(135)'],
    transitionState = 'TS203',
    kinetics = Arrhenius(A=(6.97354e-12,'cm^3/(molecule*s)'), n=0.6, Ea=(0,'kJ/mol'), T0=(1,'K'), comment="""From training reaction 18 used for C_rad/H2/Cs;H_rad
Exact match found for rate rule [C_rad/H2/Cs;H_rad]
Euclidian distance = 0
Multiplied by reaction path degeneracy 2.0
family: R_Recombination
Ea raised from -3.3 to 0 kJ/mol."""),
)

reaction(
    label = 'reaction204',
    reactants = ['[CH2]CCCC(5)', 'CH2(S)(14)'],
    products = ['[CH2]CC(C)CC(135)'],
    transitionState = 'TS204',
    kinetics = Arrhenius(A=(873476,'m^3/(mol*s)'), n=0.189, Ea=(0,'kJ/mol'), T0=(1,'K'), comment="""Estimated using template [carbene;Cs_H] for rate rule [carbene;C/H2/NonDeC]
Euclidian distance = 2.0
Multiplied by reaction path degeneracy 2.0
family: 1,2_Insertion_carbene
Ea raised from -1.5 to 0 kJ/mol."""),
)

reaction(
    label = 'reaction205',
    reactants = ['CH2(S)(14)', '[CH2]CC(C)C(120)'],
    products = ['[CH2]CC(C)CC(135)'],
    transitionState = 'TS205',
    kinetics = Arrhenius(A=(2.62043e+06,'m^3/(mol*s)'), n=0.189, Ea=(0,'kJ/mol'), T0=(1,'K'), comment="""Estimated using template [carbene;C_pri] for rate rule [carbene;C_pri/NonDeC]
Euclidian distance = 1.0
Multiplied by reaction path degeneracy 6.0
family: 1,2_Insertion_carbene
Ea raised from -1.5 to 0 kJ/mol."""),
)

reaction(
    label = 'reaction206',
    reactants = ['[CH2]C(C)CC(103)', 'CH2(T)(19)'],
    products = ['[CH2]CC(C)CC(135)'],
    transitionState = 'TS206',
    kinetics = Arrhenius(A=(759150,'m^3/(mol*s)'), n=0.575199, Ea=(34.3157,'kJ/mol'), T0=(1,'K'), comment="""Estimated using template [Y_rad;Birad] for rate rule [C_rad/H2/Cs;Birad]
Euclidian distance = 3.0
family: Birad_R_Recombination"""),
)

reaction(
    label = 'reaction207',
    reactants = ['[CH]CC(C)CC(169)', 'H(8)'],
    products = ['[CH2]CC(C)CC(135)'],
    transitionState = 'TS207',
    kinetics = Arrhenius(A=(1e+07,'m^3/(mol*s)'), n=0, Ea=(0,'kJ/mol'), T0=(1,'K'), comment="""Estimated using an average for rate rule [H_rad;Birad]
Euclidian distance = 0
family: Birad_R_Recombination"""),
)

reaction(
    label = 'reaction208',
    reactants = ['H(8)', 'C=C(C)C(C)C(170)'],
    products = ['[CH2]C(C)C(C)C(130)'],
    transitionState = 'TS208',
    kinetics = Arrhenius(A=(0.00507518,'m^3/(mol*s)'), n=2.82235, Ea=(0,'kJ/mol'), T0=(1,'K'), Tmin=(303.03,'K'), Tmax=(2000,'K'), comment="""From training reaction 102 used for Cds-CsCs_Cds-HH;HJ
Exact match found for rate rule [Cds-CsCs_Cds-HH;HJ]
Euclidian distance = 0
family: R_Addition_MultipleBond
Ea raised from -4.8 to 0 kJ/mol."""),
)

reaction(
    label = 'reaction209',
    reactants = ['C=CC(36)', 'C[CH]C(37)'],
    products = ['[CH2]C(C)C(C)C(130)'],
    transitionState = 'TS209',
    kinetics = Arrhenius(A=(818,'cm^3/(mol*s)'), n=2.41, Ea=(22.3844,'kJ/mol'), T0=(1,'K'), Tmin=(300,'K'), Tmax=(1500,'K'), comment="""Matched reaction 419 propene_2 + C3H7-2 <=> C6H13-4 in R_Addition_MultipleBond/training
This reaction matched rate rule [Cds-CsH_Cds-HH;CsJ-CsCsH]
family: R_Addition_MultipleBond"""),
)

reaction(
    label = 'reaction210',
    reactants = ['[CH3](11)', 'C=CC(C)C(124)'],
    products = ['[CH2]C(C)C(C)C(130)'],
    transitionState = 'TS210',
    kinetics = Arrhenius(A=(10000,'cm^3/(mol*s)'), n=2.41, Ea=(29.7482,'kJ/mol'), T0=(1,'K'), Tmin=(300,'K'), Tmax=(1500,'K'), comment="""From training reaction 417 used for Cds-CsH_Cds-HH;CsJ-HHH
Exact match found for rate rule [Cds-CsH_Cds-HH;CsJ-HHH]
Euclidian distance = 0
family: R_Addition_MultipleBond"""),
)

reaction(
    label = 'reaction211',
    reactants = ['[CH2]C(C)C(C)C(130)'],
    products = ['C[C](C)C(C)C(171)'],
    transitionState = 'TS211',
    kinetics = Arrhenius(A=(5.265e-07,'s^-1'), n=5.639, Ea=(102.68,'kJ/mol'), T0=(1,'K'), comment="""From training reaction 38 used for R2H_S;C_rad_out_2H;Cs_H_out_Cs2
Exact match found for rate rule [R2H_S;C_rad_out_2H;Cs_H_out_Cs2]
Euclidian distance = 0
family: intra_H_migration"""),
)

reaction(
    label = 'reaction212',
    reactants = ['[CH2][CH]C(38)', 'C[CH]C(37)'],
    products = ['[CH2]C(C)C(C)C(130)'],
    transitionState = 'TS212',
    kinetics = Arrhenius(A=(2.79546e+07,'m^3/(mol*s)'), n=-0.206, Ea=(0,'kJ/mol'), T0=(1,'K'), comment="""Estimated using an average for rate rule [Y_rad;C_rad/H/NonDeC]
Euclidian distance = 0
family: R_Recombination
Ea raised from -0.0 to 0 kJ/mol."""),
)

reaction(
    label = 'reaction213',
    reactants = ['[CH3](11)', '[CH2]C(C)[CH]C(154)'],
    products = ['[CH2]C(C)C(C)C(130)'],
    transitionState = 'TS213',
    kinetics = Arrhenius(A=(6.64e+14,'cm^3/(mol*s)'), n=-0.57, Ea=(0,'kJ/mol'), T0=(1,'K'), Tmin=(713,'K'), Tmax=(1800,'K'), comment="""From training reaction 67 used for C_methyl;C_rad/H/NonDeC
Exact match found for rate rule [C_rad/H/NonDeC;C_methyl]
Euclidian distance = 0
family: R_Recombination"""),
)

reaction(
    label = 'reaction214',
    reactants = ['H(8)', '[CH2]C(C)[C](C)C(172)'],
    products = ['[CH2]C(C)C(C)C(130)'],
    transitionState = 'TS214',
    kinetics = Arrhenius(A=(6.68468e+06,'m^3/(mol*s)'), n=-0.0135, Ea=(0,'kJ/mol'), T0=(1,'K'), comment="""Estimated using template [C_rad/Cs3;Y_rad] for rate rule [C_rad/Cs3;H_rad]
Euclidian distance = 1.0
family: R_Recombination
Ea raised from -0.9 to 0 kJ/mol."""),
)

reaction(
    label = 'reaction215',
    reactants = ['[CH3](11)', '[CH2][CH]C(C)C(126)'],
    products = ['[CH2]C(C)C(C)C(130)'],
    transitionState = 'TS215',
    kinetics = Arrhenius(A=(1.66881e+08,'m^3/(mol*s)'), n=-0.401267, Ea=(0,'kJ/mol'), T0=(1,'K'), comment="""Estimated using an average for rate rule [Y_rad;C_methyl]
Euclidian distance = 0
family: R_Recombination
Ea raised from -6.7 to 0 kJ/mol."""),
)

reaction(
    label = 'reaction216',
    reactants = ['H(8)', '[CH2][C](C)C(C)C(173)'],
    products = ['[CH2]C(C)C(C)C(130)'],
    transitionState = 'TS216',
    kinetics = Arrhenius(A=(4.34078e+06,'m^3/(mol*s)'), n=0.278577, Ea=(0,'kJ/mol'), T0=(1,'K'), comment="""Estimated using an average for rate rule [Y_rad;H_rad]
Euclidian distance = 0
family: R_Recombination
Ea raised from -1.4 to 0 kJ/mol."""),
)

reaction(
    label = 'reaction217',
    reactants = ['H(8)', '[CH2]C(C)C([CH2])C(174)'],
    products = ['[CH2]C(C)C(C)C(130)'],
    transitionState = 'TS217',
    kinetics = Arrhenius(A=(6.97354e-12,'cm^3/(molecule*s)'), n=0.6, Ea=(0,'kJ/mol'), T0=(1,'K'), comment="""From training reaction 18 used for C_rad/H2/Cs;H_rad
Exact match found for rate rule [C_rad/H2/Cs;H_rad]
Euclidian distance = 0
Multiplied by reaction path degeneracy 2.0
family: R_Recombination
Ea raised from -3.3 to 0 kJ/mol."""),
)

reaction(
    label = 'reaction218',
    reactants = ['H(8)', '[CH2]C([CH2])C(C)C(175)'],
    products = ['[CH2]C(C)C(C)C(130)'],
    transitionState = 'TS218',
    kinetics = Arrhenius(A=(6.97354e-12,'cm^3/(molecule*s)'), n=0.6, Ea=(0,'kJ/mol'), T0=(1,'K'), comment="""From training reaction 18 used for C_rad/H2/Cs;H_rad
Exact match found for rate rule [C_rad/H2/Cs;H_rad]
Euclidian distance = 0
Multiplied by reaction path degeneracy 2.0
family: R_Recombination
Ea raised from -3.3 to 0 kJ/mol."""),
)

reaction(
    label = 'reaction219',
    reactants = ['[CH2]C(C)CC(103)', 'CH2(S)(14)'],
    products = ['[CH2]C(C)C(C)C(130)'],
    transitionState = 'TS219',
    kinetics = Arrhenius(A=(873476,'m^3/(mol*s)'), n=0.189, Ea=(0,'kJ/mol'), T0=(1,'K'), comment="""Estimated using template [carbene;Cs_H] for rate rule [carbene;C/H2/NonDeC]
Euclidian distance = 2.0
Multiplied by reaction path degeneracy 2.0
family: 1,2_Insertion_carbene
Ea raised from -1.5 to 0 kJ/mol."""),
)

reaction(
    label = 'reaction220',
    reactants = ['CH2(S)(14)', '[CH2]CC(C)C(120)'],
    products = ['[CH2]C(C)C(C)C(130)'],
    transitionState = 'TS220',
    kinetics = Arrhenius(A=(143764,'m^3/(mol*s)'), n=0.444, Ea=(0,'kJ/mol'), T0=(1,'K'), comment="""Estimated using an average for rate rule [carbene;R_H]
Euclidian distance = 0
Multiplied by reaction path degeneracy 2.0
family: 1,2_Insertion_carbene
Ea raised from -5.1 to 0 kJ/mol."""),
)

reaction(
    label = 'reaction221',
    reactants = ['CH2(T)(19)', 'C[CH]C(C)C(128)'],
    products = ['[CH2]C(C)C(C)C(130)'],
    transitionState = 'TS221',
    kinetics = Arrhenius(A=(759150,'m^3/(mol*s)'), n=0.575199, Ea=(34.3157,'kJ/mol'), T0=(1,'K'), comment="""Estimated using template [Y_rad;Birad] for rate rule [C_rad/H/NonDeC;Birad]
Euclidian distance = 3.0
family: Birad_R_Recombination"""),
)

reaction(
    label = 'reaction222',
    reactants = ['H(8)', '[CH]C(C)C(C)C(176)'],
    products = ['[CH2]C(C)C(C)C(130)'],
    transitionState = 'TS222',
    kinetics = Arrhenius(A=(1e+07,'m^3/(mol*s)'), n=0, Ea=(0,'kJ/mol'), T0=(1,'K'), comment="""Estimated using an average for rate rule [H_rad;Birad]
Euclidian distance = 0
family: Birad_R_Recombination"""),
)

reaction(
    label = 'reaction223',
    reactants = ['H(8)', 'C=CC(C)(C)C(177)'],
    products = ['C[CH]C(C)(C)C(160)'],
    transitionState = 'TS223',
    kinetics = Arrhenius(A=(3.36e+08,'cm^3/(mol*s)'), n=1.56, Ea=(2.5104,'kJ/mol'), T0=(1,'K'), comment="""From training reaction 9 used for Cds-HH_Cds-CsH;HJ
Exact match found for rate rule [Cds-HH_Cds-CsH;HJ]
Euclidian distance = 0
family: R_Addition_MultipleBond"""),
)

reaction(
    label = 'reaction224',
    reactants = ['[CH3](11)', 'CC=C(C)C(178)'],
    products = ['C[CH]C(C)(C)C(160)'],
    transitionState = 'TS224',
    kinetics = Arrhenius(A=(6260,'cm^3/(mol*s)'), n=2.41, Ea=(32.3423,'kJ/mol'), T0=(1,'K'), Tmin=(300,'K'), Tmax=(1500,'K'), comment="""Matched reaction 633 CH3 + C5H10-2 <=> C6H13-8 in R_Addition_MultipleBond/training
This reaction matched rate rule [Cds-CsCs_Cds-CsH;CsJ-HHH]
family: R_Addition_MultipleBond"""),
)

reaction(
    label = 'reaction225',
    reactants = ['C[CH]C(C)(C)C(160)'],
    products = ['[CH2]CC(C)(C)C(161)'],
    transitionState = 'TS225',
    kinetics = Arrhenius(A=(5.04e-11,'s^-1'), n=6.833, Ea=(117.248,'kJ/mol'), T0=(1,'K'), comment="""From training reaction 33 used for R2H_S;C_rad_out_H/NonDeC;Cs_H_out_2H
Exact match found for rate rule [R2H_S;C_rad_out_H/NonDeC;Cs_H_out_2H]
Euclidian distance = 0
Multiplied by reaction path degeneracy 3.0
family: intra_H_migration"""),
)

reaction(
    label = 'reaction226',
    reactants = ['[CH3](11)', 'C[CH][C](C)C(179)'],
    products = ['C[CH]C(C)(C)C(160)'],
    transitionState = 'TS226',
    kinetics = Arrhenius(A=(1.66881e+08,'m^3/(mol*s)'), n=-0.401267, Ea=(0,'kJ/mol'), T0=(1,'K'), comment="""Estimated using an average for rate rule [Y_rad;C_methyl]
Euclidian distance = 0
family: R_Recombination
Ea raised from -6.7 to 0 kJ/mol."""),
)

reaction(
    label = 'reaction227',
    reactants = ['H(8)', '[CH2]C(C)(C)[CH]C(163)'],
    products = ['C[CH]C(C)(C)C(160)'],
    transitionState = 'TS227',
    kinetics = Arrhenius(A=(3.48677e-12,'cm^3/(molecule*s)'), n=0.6, Ea=(0,'kJ/mol'), T0=(1,'K'), comment="""From training reaction 18 used for C_rad/H2/Cs;H_rad
Exact match found for rate rule [C_rad/H2/Cs;H_rad]
Euclidian distance = 0
family: R_Recombination
Ea raised from -3.3 to 0 kJ/mol."""),
)

reaction(
    label = 'reaction228',
    reactants = ['[CH2][CH]C(C)(C)C(180)', 'H(8)'],
    products = ['C[CH]C(C)(C)C(160)'],
    transitionState = 'TS228',
    kinetics = Arrhenius(A=(4.34078e+06,'m^3/(mol*s)'), n=0.278577, Ea=(0,'kJ/mol'), T0=(1,'K'), comment="""Estimated using an average for rate rule [Y_rad;H_rad]
Euclidian distance = 0
family: R_Recombination
Ea raised from -1.4 to 0 kJ/mol."""),
)

reaction(
    label = 'reaction229',
    reactants = ['CH2(S)(14)', 'C[CH]C(C)C(128)'],
    products = ['C[CH]C(C)(C)C(160)'],
    transitionState = 'TS229',
    kinetics = Arrhenius(A=(71881.9,'m^3/(mol*s)'), n=0.444, Ea=(0,'kJ/mol'), T0=(1,'K'), comment="""Estimated using an average for rate rule [carbene;R_H]
Euclidian distance = 0
family: 1,2_Insertion_carbene
Ea raised from -5.1 to 0 kJ/mol."""),
)

reaction(
    label = 'reaction230',
    reactants = ['C[CH]C(C)(C)C(160)'],
    products = ['C[C](C)C(C)C(171)'],
    transitionState = 'TS230',
    kinetics = Arrhenius(A=(3.99e+08,'s^-1'), n=1.36, Ea=(157.318,'kJ/mol'), T0=(1,'K'), comment="""Estimated using template [cCs(-R!HR!H)CJ;CsJ;CH3] for rate rule [cCs(-R!HR!H)CJ;CsJ-CsH;CH3]
Euclidian distance = 1.0
Multiplied by reaction path degeneracy 3.0
family: 1,2_shiftC"""),
)

reaction(
    label = 'reaction231',
    reactants = ['C[C](C)C(181)', '[CH]C(20)'],
    products = ['C[CH]C(C)(C)C(160)'],
    transitionState = 'TS231',
    kinetics = Arrhenius(A=(759150,'m^3/(mol*s)'), n=0.575199, Ea=(34.3157,'kJ/mol'), T0=(1,'K'), comment="""Estimated using template [Y_rad;Birad] for rate rule [C_rad/Cs3;Birad]
Euclidian distance = 4.0
family: Birad_R_Recombination"""),
)

reaction(
    label = 'reaction232',
    reactants = ['[CH3](11)', '[CH]C(C)(C)C(182)'],
    products = ['C[CH]C(C)(C)C(160)'],
    transitionState = 'TS232',
    kinetics = Arrhenius(A=(759150,'m^3/(mol*s)'), n=0.575199, Ea=(34.3157,'kJ/mol'), T0=(1,'K'), comment="""Estimated using template [Y_rad;Birad] for rate rule [C_methyl;Birad]
Euclidian distance = 2.0
family: Birad_R_Recombination"""),
)

reaction(
    label = 'reaction233',
    reactants = ['H(8)', 'C[C]C(C)(C)C(183)'],
    products = ['C[CH]C(C)(C)C(160)'],
    transitionState = 'TS233',
    kinetics = Arrhenius(A=(1e+07,'m^3/(mol*s)'), n=0, Ea=(0,'kJ/mol'), T0=(1,'K'), comment="""Estimated using an average for rate rule [H_rad;Birad]
Euclidian distance = 0
family: Birad_R_Recombination"""),
)

reaction(
    label = 'reaction234',
    reactants = ['H(8)', 'C=CC(C)(C)C(177)'],
    products = ['[CH2]CC(C)(C)C(161)'],
    transitionState = 'TS234',
    kinetics = Arrhenius(A=(3.01e+08,'cm^3/(mol*s)'), n=1.6, Ea=(10.0416,'kJ/mol'), T0=(1,'K'), comment="""From training reaction 10 used for Cds-CsH_Cds-HH;HJ
Exact match found for rate rule [Cds-CsH_Cds-HH;HJ]
Euclidian distance = 0
family: R_Addition_MultipleBond"""),
)

reaction(
    label = 'reaction235',
    reactants = ['C[C](C)C(181)', 'C=C(17)'],
    products = ['[CH2]CC(C)(C)C(161)'],
    transitionState = 'TS235',
    kinetics = Arrhenius(A=(2520,'cm^3/(mol*s)'), n=2.41, Ea=(10.2508,'kJ/mol'), T0=(1,'K'), Tmin=(300,'K'), Tmax=(1500,'K'), comment="""Matched reaction 222 C2H4 + C4H9 <=> C6H13 in R_Addition_MultipleBond/training
This reaction matched rate rule [Cds-HH_Cds-HH;CsJ-CsCsCs]
family: R_Addition_MultipleBond"""),
)

reaction(
    label = 'reaction236',
    reactants = ['C[C](C)C(181)', '[CH2][CH2](18)'],
    products = ['[CH2]CC(C)(C)C(161)'],
    transitionState = 'TS236',
    kinetics = Arrhenius(A=(1.33694e+07,'m^3/(mol*s)'), n=-0.0135, Ea=(0,'kJ/mol'), T0=(1,'K'), comment="""Estimated using an average for rate rule [Y_rad;C_rad/Cs3]
Euclidian distance = 0
Multiplied by reaction path degeneracy 2.0
family: R_Recombination
Ea raised from -0.9 to 0 kJ/mol."""),
)

reaction(
    label = 'reaction237',
    reactants = ['[CH2]C[C](C)C(108)', '[CH3](11)'],
    products = ['[CH2]CC(C)(C)C(161)'],
    transitionState = 'TS237',
    kinetics = Arrhenius(A=(4.88e+15,'cm^3/(mol*s)','*|/',2), n=-1, Ea=(0,'kJ/mol'), T0=(1,'K'), Tmin=(300,'K'), Tmax=(2500,'K'), comment="""From training reaction 68 used for C_methyl;C_rad/Cs3
Exact match found for rate rule [C_rad/Cs3;C_methyl]
Euclidian distance = 0
family: R_Recombination"""),
)

reaction(
    label = 'reaction238',
    reactants = ['[CH2][CH]C(C)(C)C(180)', 'H(8)'],
    products = ['[CH2]CC(C)(C)C(161)'],
    transitionState = 'TS238',
    kinetics = Arrhenius(A=(4.34078e+06,'m^3/(mol*s)'), n=0.278577, Ea=(0,'kJ/mol'), T0=(1,'K'), comment="""Estimated using an average for rate rule [Y_rad;H_rad]
Euclidian distance = 0
family: R_Recombination
Ea raised from -1.4 to 0 kJ/mol."""),
)

reaction(
    label = 'reaction239',
    reactants = ['H(8)', '[CH2]CC([CH2])(C)C(165)'],
    products = ['[CH2]CC(C)(C)C(161)'],
    transitionState = 'TS239',
    kinetics = Arrhenius(A=(3.48677e-12,'cm^3/(molecule*s)'), n=0.6, Ea=(0,'kJ/mol'), T0=(1,'K'), comment="""From training reaction 18 used for C_rad/H2/Cs;H_rad
Exact match found for rate rule [C_rad/H2/Cs;H_rad]
Euclidian distance = 0
family: R_Recombination
Ea raised from -3.3 to 0 kJ/mol."""),
)

reaction(
    label = 'reaction240',
    reactants = ['CH2(S)(14)', '[CH2]CC(C)C(120)'],
    products = ['[CH2]CC(C)(C)C(161)'],
    transitionState = 'TS240',
    kinetics = Arrhenius(A=(436738,'m^3/(mol*s)'), n=0.189, Ea=(0,'kJ/mol'), T0=(1,'K'), comment="""Estimated using template [carbene;Cs_H] for rate rule [carbene;C/H/Cs3]
Euclidian distance = 3.0
family: 1,2_Insertion_carbene
Ea raised from -1.5 to 0 kJ/mol."""),
)

reaction(
    label = 'reaction241',
    reactants = ['CH2(T)(19)', '[CH2]C(C)(C)C(166)'],
    products = ['[CH2]CC(C)(C)C(161)'],
    transitionState = 'TS241',
    kinetics = Arrhenius(A=(759150,'m^3/(mol*s)'), n=0.575199, Ea=(34.3157,'kJ/mol'), T0=(1,'K'), comment="""Estimated using template [Y_rad;Birad] for rate rule [C_rad/H2/Cs;Birad]
Euclidian distance = 3.0
family: Birad_R_Recombination"""),
)

reaction(
    label = 'reaction242',
    reactants = ['H(8)', '[CH]CC(C)(C)C(184)'],
    products = ['[CH2]CC(C)(C)C(161)'],
    transitionState = 'TS242',
    kinetics = Arrhenius(A=(1e+07,'m^3/(mol*s)'), n=0, Ea=(0,'kJ/mol'), T0=(1,'K'), comment="""Estimated using an average for rate rule [H_rad;Birad]
Euclidian distance = 0
family: Birad_R_Recombination"""),
)

network(
    label = '50',
    isomers = [
        'C[CH]CCCC(49)',
        '[CH2]CCCCC(10)',
        'CC[CH]CCC(50)',
        '[CH2]C(C)CCC(68)',
        'CCC[C](C)C(92)',
        '[CH2]CCC(C)C(95)',
        'CC[CH]C(C)C(93)',
        '[CH2]C(CC)CC(83)',
        'C[CH]CC(C)C(94)',
        'CC[C](C)CC(134)',
        'C[CH]C(C)CC(129)',
        '[CH2]C(C)(C)CC(113)',
        '[CH2]CC(C)CC(135)',
        '[CH2]C(C)C(C)C(130)',
        'C[CH]C(C)(C)C(160)',
        '[CH2]CC(C)(C)C(161)',
    ],
    reactants = [
        ('[CH2]C[CH]C(44)', 'C[CH2](6)'),
        ('[CH2]CCC(3)', 'C=C(17)'),
    ],
    bathGas = {
        'Ne': 0.25,
        'Ar': 0.25,
        'He': 0.25,
        'N2': 0.25,
    },
)

pressureDependence(
    label = '50',
    Tmin = (1200,'K'),
    Tmax = (1500,'K'),
    Tcount = 10,
    Tlist = ([1201.48,1213.22,1236.21,1269.31,1310.55,1356.92,1404.16,1447.02,1479.84,1497.7],'K'),
    Pmin = (1,'atm'),
    Pmax = (10,'atm'),
    Pcount = 10,
    Plist = ([1.02771,1.14872,1.41959,1.89986,2.67608,3.83649,5.40396,7.23219,8.93758,9.98989],'bar'),
    maximumGrainSize = (0.5,'kcal/mol'),
    minimumGrainCount = 250,
    method = 'modified strong collision',
    interpolationModel = ('Chebyshev', 6, 4),
    activeKRotor = True,
    activeJRotor = True,
    rmgmode = True,
)

