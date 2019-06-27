species(
    label = '[CH2][CH]CC([O])O[O](8173)',
    structure = SMILES('[CH2][CH]CC([O])O[O]'),
    E0 = (344.545,'kJ/mol'),
    modes = [
        HarmonicOscillator(frequencies=([3025,407.5,1350,352.5,492.5,1135,1000,3000,3100,440,815,1455,1000,1380,1390,370,380,2900,435,2750,2850,1437.5,1250,1305,750,350,180,1676.47,1676.5],'cm^-1')),
        HinderedRotor(inertia=(0.00236438,'amu*angstrom^2'), symmetry=1, barrier=(4.71718,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.204963,'amu*angstrom^2'), symmetry=1, barrier=(4.7125,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.205018,'amu*angstrom^2'), symmetry=1, barrier=(4.71378,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.205013,'amu*angstrom^2'), symmetry=1, barrier=(4.71365,'kJ/mol'), semiclassical=False),
    ],
    spinMultiplicity = 5,
    opticalIsomers = 1,
    molecularWeight = (102.089,'amu'),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[1.24374,0.072974,-0.00012519,1.25511e-07,-4.76779e-11,41526.5,32.2605], Tmin=(100,'K'), Tmax=(851.686,'K')), NASAPolynomial(coeffs=[1.47381,0.042879,-2.10857e-05,4.02209e-09,-2.75423e-13,42539.6,37.3653], Tmin=(851.686,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(344.545,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(291.007,'J/(mol*K)'), comment="""Thermo group additivity estimation: group(O2s-OsCs) + group(O2s-CsH) + group(O2s-OsH) + group(Cs-CsCsHH) + group(Cs-CsCsHH) + group(Cs-CsOsOsH) + group(Cs-CsHHH) + radical(ROOJ) + radical(RCCJ) + radical(CCOJ) + radical(RCCJC)"""),
)

species(
    label = '[CH2]C=C(87)',
    structure = SMILES('[CH2]C=C'),
    E0 = (157.623,'kJ/mol'),
    modes = [
        HarmonicOscillator(frequencies=([2950,3100,1380,975,1025,1650,3010,987.5,1337.5,450,1655,3000,3100,440,815,1455,1000],'cm^-1')),
        HinderedRotor(inertia=(0.570287,'amu*angstrom^2'), symmetry=1, barrier=(32.8573,'kJ/mol'), semiclassical=False),
    ],
    spinMultiplicity = 2,
    opticalIsomers = 1,
    molecularWeight = (41.0718,'amu'),
    collisionModel = TransportData(shapeIndex=2, epsilon=(2161.77,'J/mol'), sigma=(4.85,'angstroms'), dipoleMoment=(0,'C*m'), polarizability=(0,'angstroms^3'), rotrelaxcollnum=1.0, comment="""GRI-Mech"""),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[3.3193,0.00566487,4.27449e-05,-5.78831e-08,2.21699e-11,18990.6,9.19646], Tmin=(100,'K'), Tmax=(951.999,'K')), NASAPolynomial(coeffs=[7.55715,0.0114811,-3.63952e-06,6.63584e-10,-4.95318e-14,17113.3,-16.6624], Tmin=(951.999,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(157.623,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(178.761,'J/(mol*K)'), comment="""Thermo group additivity estimation: group(Cs-(Cds-Cds)HHH) + group(Cds-CdsCsH) + group(Cds-CdsHH) + radical(Allyl_P)"""),
)

species(
    label = '[O]OC=O(5472)',
    structure = SMILES('[O]OC=O'),
    E0 = (-195.495,'kJ/mol'),
    modes = [
        HarmonicOscillator(frequencies=([492.5,1135,1000,739.225,739.248,739.254,739.261,739.261],'cm^-1')),
        HinderedRotor(inertia=(0.00030847,'amu*angstrom^2'), symmetry=1, barrier=(0.119627,'kJ/mol'), semiclassical=False),
    ],
    spinMultiplicity = 2,
    opticalIsomers = 1,
    molecularWeight = (61.0168,'amu'),
    collisionModel = TransportData(shapeIndex=2, epsilon=(3570.08,'J/mol'), sigma=(5.61676,'angstroms'), dipoleMoment=(0,'C*m'), polarizability=(0,'angstroms^3'), rotrelaxcollnum=0, comment="""Epsilon & sigma estimated with Tc=557.64 K, Pc=45.72 bar (from Joback method)"""),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[3.38375,0.00191416,4.59985e-05,-6.61808e-08,2.67351e-11,-23479.8,12.2589], Tmin=(100,'K'), Tmax=(935.456,'K')), NASAPolynomial(coeffs=[10.7708,0.000103189,1.15693e-06,-1.97269e-10,7.46894e-15,-26164.6,-29.8498], Tmin=(935.456,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(-195.495,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(103.931,'J/(mol*K)'), comment="""Thermo group additivity estimation: group(O2s-O2s(Cds-O2d)) + group(O2s-OsH) + group(Cds-OdOsH) + radical(C(=O)OOJ)"""),
)

species(
    label = 'H(8)',
    structure = SMILES('[H]'),
    E0 = (211.805,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
    molecularWeight = (1.00794,'amu'),
    collisionModel = TransportData(shapeIndex=0, epsilon=(1205.6,'J/mol'), sigma=(2.05,'angstroms'), dipoleMoment=(0,'C*m'), polarizability=(0,'angstroms^3'), rotrelaxcollnum=0.0, comment="""GRI-Mech"""),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[2.5,9.24385e-15,-1.3678e-17,6.66185e-21,-1.00107e-24,25474.2,-0.444973], Tmin=(100,'K'), Tmax=(3459.6,'K')), NASAPolynomial(coeffs=[2.5,9.20456e-12,-3.58608e-15,6.15199e-19,-3.92042e-23,25474.2,-0.444973], Tmin=(3459.6,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(211.805,'kJ/mol'), Cp0=(20.7862,'J/(mol*K)'), CpInf=(20.7862,'J/(mol*K)'), label="""H""", comment="""Thermo library: primaryThermoLibrary"""),
)

species(
    label = 'C=C[CH]C([O])O[O](8518)',
    structure = SMILES('C=C[CH]C([O])O[O]'),
    E0 = (189.551,'kJ/mol'),
    modes = [
        HarmonicOscillator(frequencies=([492.5,1135,1000,2950,3100,1380,975,1025,1650,3010,987.5,1337.5,450,1655,1380,1390,370,380,2900,435,3025,407.5,1350,352.5,247.329,247.339,1680.38],'cm^-1')),
        HinderedRotor(inertia=(0.270626,'amu*angstrom^2'), symmetry=1, barrier=(11.5076,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(1.01617,'amu*angstrom^2'), symmetry=1, barrier=(43.3833,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(1.01248,'amu*angstrom^2'), symmetry=1, barrier=(43.383,'kJ/mol'), semiclassical=False),
    ],
    spinMultiplicity = 4,
    opticalIsomers = 1,
    molecularWeight = (101.081,'amu'),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[1.46346,0.0629289,-6.9289e-05,3.59254e-08,-1.97701e-12,22882.3,24.6775], Tmin=(100,'K'), Tmax=(572.353,'K')), NASAPolynomial(coeffs=[7.02067,0.0334739,-1.66843e-05,3.29394e-09,-2.34229e-13,22092.5,-0.374396], Tmin=(572.353,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(189.551,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(270.22,'J/(mol*K)'), comment="""Thermo group additivity estimation: group(O2s-OsCs) + group(O2s-CsH) + group(O2s-OsH) + group(Cs-(Cds-Cds)CsHH) + group(Cs-CsOsOsH) + group(Cds-CdsCsH) + group(Cds-CdsHH) + radical(C=CCJCO) + radical(CCOJ) + radical(ROOJ)"""),
)

species(
    label = 'C=CC[C]([O])O[O](5392)',
    structure = SMILES('C=CC[C]([O])O[O]'),
    E0 = (277.881,'kJ/mol'),
    modes = [
        HarmonicOscillator(frequencies=([492.5,1135,1000,360,370,350,2950,3100,1380,975,1025,1650,3010,987.5,1337.5,450,1655,2750,2850,1437.5,1250,1305,750,350,180,180,1611.29],'cm^-1')),
        HinderedRotor(inertia=(0.200123,'amu*angstrom^2'), symmetry=1, barrier=(4.60121,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.195767,'amu*angstrom^2'), symmetry=1, barrier=(4.50108,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.197702,'amu*angstrom^2'), symmetry=1, barrier=(4.54555,'kJ/mol'), semiclassical=False),
    ],
    spinMultiplicity = 4,
    opticalIsomers = 1,
    molecularWeight = (101.081,'amu'),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[1.35373,0.0649565,-9.60355e-05,8.72158e-08,-3.21097e-11,33510.1,28.5067], Tmin=(100,'K'), Tmax=(797.789,'K')), NASAPolynomial(coeffs=[5.23015,0.0341917,-1.68909e-05,3.27926e-09,-2.28989e-13,33252.1,12.9405], Tmin=(797.789,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(277.881,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(270.22,'J/(mol*K)'), comment="""Thermo group additivity estimation: group(O2s-OsCs) + group(O2s-CsH) + group(O2s-OsH) + group(Cs-(Cds-Cds)CsHH) + group(Cs-CsOsOsH) + group(Cds-CdsCsH) + group(Cds-CdsHH) + radical(CCOJ) + radical(Cs_P) + radical(ROOJ)"""),
)

species(
    label = '[O][CH]O[O](8201)',
    structure = SMILES('[O][CH]O[O]'),
    E0 = (228.888,'kJ/mol'),
    modes = [
        HarmonicOscillator(frequencies=([492.5,1135,1000,3025,407.5,1350,352.5,1824.75],'cm^-1')),
        HinderedRotor(inertia=(0.270955,'amu*angstrom^2'), symmetry=1, barrier=(6.2298,'kJ/mol'), semiclassical=False),
    ],
    spinMultiplicity = 4,
    opticalIsomers = 1,
    molecularWeight = (61.0168,'amu'),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[2.97073,0.0345271,-8.63149e-05,9.83474e-08,-3.87219e-11,27554.6,13.3771], Tmin=(100,'K'), Tmax=(878.005,'K')), NASAPolynomial(coeffs=[-0.738633,0.0210949,-1.15487e-05,2.23205e-09,-1.51294e-13,29375,37.4477], Tmin=(878.005,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(228.888,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(103.931,'J/(mol*K)'), comment="""Thermo group additivity estimation: group(O2s-OsCs) + group(O2s-CsH) + group(O2s-OsH) + group(Cs-OsOsHH) + radical(OCOJ) + radical(ROOJ) + radical(OCJO)"""),
)

species(
    label = '[CH2][CH][CH2](6136)',
    structure = SMILES('[CH2][CH][CH2]'),
    E0 = (484.293,'kJ/mol'),
    modes = [
        HarmonicOscillator(frequencies=([3025,407.5,1350,352.5,3000,3033.33,3066.67,3100,415,465,780,850,1435,1475,900,1100],'cm^-1')),
        HinderedRotor(inertia=(0.00132719,'amu*angstrom^2'), symmetry=1, barrier=(2.41051,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.00132749,'amu*angstrom^2'), symmetry=1, barrier=(2.41088,'kJ/mol'), semiclassical=False),
    ],
    spinMultiplicity = 4,
    opticalIsomers = 1,
    molecularWeight = (41.0718,'amu'),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[3.34321,0.013802,2.16426e-06,-5.76329e-09,1.61332e-12,58271,14.955], Tmin=(100,'K'), Tmax=(1447.11,'K')), NASAPolynomial(coeffs=[4.39505,0.0167645,-6.99091e-06,1.25741e-09,-8.38108e-14,57351.9,7.36811], Tmin=(1447.11,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(484.293,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(174.604,'J/(mol*K)'), comment="""Thermo group additivity estimation: group(Cs-CsCsHH) + group(Cs-CsHHH) + group(Cs-CsHHH) + radical(RCCJ) + radical(CCJC) + radical(RCCJ)"""),
)

species(
    label = 'O2(2)',
    structure = SMILES('[O][O]'),
    E0 = (-8.62683,'kJ/mol'),
    modes = [
        HarmonicOscillator(frequencies=([1487.4],'cm^-1')),
    ],
    spinMultiplicity = 3,
    opticalIsomers = 1,
    molecularWeight = (31.9988,'amu'),
    collisionModel = TransportData(shapeIndex=1, epsilon=(892.977,'J/mol'), sigma=(3.458,'angstroms'), dipoleMoment=(0,'C*m'), polarizability=(1.6,'angstroms^3'), rotrelaxcollnum=3.8, comment="""GRI-Mech"""),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[3.53732,-0.00121571,5.31618e-06,-4.89443e-09,1.45845e-12,-1038.59,4.68368], Tmin=(100,'K'), Tmax=(1074.56,'K')), NASAPolynomial(coeffs=[3.15382,0.00167804,-7.69971e-07,1.51275e-10,-1.08782e-14,-1040.82,6.16754], Tmin=(1074.56,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(-8.62683,'kJ/mol'), Cp0=(29.1007,'J/(mol*K)'), CpInf=(37.4151,'J/(mol*K)'), label="""O2""", comment="""Thermo library: primaryThermoLibrary"""),
)

species(
    label = 'C=CC[CH][O](692)',
    structure = SMILES('C=CC[CH][O]'),
    E0 = (229.654,'kJ/mol'),
    modes = [
        HarmonicOscillator(frequencies=([2750,2850,1437.5,1250,1305,750,350,2950,3100,1380,975,1025,1650,3010,987.5,1337.5,450,1655,3025,407.5,1350,352.5,245.72,245.721,1409.84],'cm^-1')),
        HinderedRotor(inertia=(0.128783,'amu*angstrom^2'), symmetry=1, barrier=(5.51785,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.128782,'amu*angstrom^2'), symmetry=1, barrier=(5.51783,'kJ/mol'), semiclassical=False),
    ],
    spinMultiplicity = 3,
    opticalIsomers = 1,
    molecularWeight = (70.0898,'amu'),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[2.16347,0.0442974,-5.21503e-05,4.48062e-08,-1.67944e-11,27683.3,19.7842], Tmin=(100,'K'), Tmax=(774.394,'K')), NASAPolynomial(coeffs=[3.80669,0.0302192,-1.40523e-05,2.68593e-09,-1.87067e-13,27596.5,13.359], Tmin=(774.394,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(229.654,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(249.434,'J/(mol*K)'), comment="""Thermo group additivity estimation: group(O2s-CsH) + group(Cs-(Cds-Cds)CsHH) + group(Cs-CsOsHH) + group(Cds-CdsCsH) + group(Cds-CdsHH) + radical(CCsJOH) + radical(CCOJ)"""),
)

species(
    label = '[CH2]C[CH]C([O])O[O](8519)',
    structure = SMILES('[CH2]C[CH]C([O])O[O]'),
    E0 = (350.514,'kJ/mol'),
    modes = [
        HarmonicOscillator(frequencies=([3025,407.5,1350,352.5,492.5,1135,1000,3000,3100,440,815,1455,1000,1380,1390,370,380,2900,435,2750,2850,1437.5,1250,1305,750,350,180,1615.58,1616.53],'cm^-1')),
        HinderedRotor(inertia=(0.167635,'amu*angstrom^2'), symmetry=1, barrier=(3.85425,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.167171,'amu*angstrom^2'), symmetry=1, barrier=(3.84359,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.167562,'amu*angstrom^2'), symmetry=1, barrier=(3.85259,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.167216,'amu*angstrom^2'), symmetry=1, barrier=(3.84462,'kJ/mol'), semiclassical=False),
    ],
    spinMultiplicity = 5,
    opticalIsomers = 1,
    molecularWeight = (102.089,'amu'),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[1.04891,0.0750437,-0.000122529,1.17259e-07,-4.36464e-11,42253.5,32.5293], Tmin=(100,'K'), Tmax=(833.516,'K')), NASAPolynomial(coeffs=[3.9408,0.0395113,-1.96152e-05,3.7776e-09,-2.60988e-13,42523.6,23.6165], Tmin=(833.516,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(350.514,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(291.007,'J/(mol*K)'), comment="""Thermo group additivity estimation: group(O2s-OsCs) + group(O2s-CsH) + group(O2s-OsH) + group(Cs-CsCsHH) + group(Cs-CsCsHH) + group(Cs-CsOsOsH) + group(Cs-CsHHH) + radical(ROOJ) + radical(CCOJ) + radical(CCJCOOH) + radical(RCCJ)"""),
)

species(
    label = '[CH2][CH]C[C](O)O[O](8520)',
    structure = SMILES('[CH2][CH]C[C](O)O[O]'),
    E0 = (324.086,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
    molecularWeight = (102.089,'amu'),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[0.933853,0.0762527,-0.000122508,1.11333e-07,-3.92892e-11,39080.5,33.2015], Tmin=(100,'K'), Tmax=(854.904,'K')), NASAPolynomial(coeffs=[5.97638,0.0341374,-1.61146e-05,3.02382e-09,-2.05318e-13,38895.2,13.6234], Tmin=(854.904,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(324.086,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(286.849,'J/(mol*K)'), comment="""Thermo group additivity estimation: group(O2s-OsCs) + group(O2s-CsH) + group(O2s-OsH) + group(Cs-CsCsHH) + group(Cs-CsCsHH) + group(Cs-CsOsOsH) + group(Cs-CsHHH) + radical(Cs_P) + radical(RCCJC) + radical(ROOJ) + radical(RCCJ)"""),
)

species(
    label = '[CH2][CH]C[C]([O])OO(1366)',
    structure = SMILES('[CH2][CH]C[C]([O])OO'),
    E0 = (397.787,'kJ/mol'),
    modes = [
        HarmonicOscillator(frequencies=([3025,407.5,1350,352.5,3615,1310,387.5,850,1000,360,370,350,3000,3100,440,815,1455,1000,2750,2850,1437.5,1250,1305,750,350,260.959,915.613,2087.67],'cm^-1')),
        HinderedRotor(inertia=(0.086067,'amu*angstrom^2'), symmetry=1, barrier=(3.36896,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.086067,'amu*angstrom^2'), symmetry=1, barrier=(3.36896,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.086067,'amu*angstrom^2'), symmetry=1, barrier=(3.36896,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.086067,'amu*angstrom^2'), symmetry=1, barrier=(3.36896,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.086067,'amu*angstrom^2'), symmetry=1, barrier=(3.36896,'kJ/mol'), semiclassical=False),
    ],
    spinMultiplicity = 5,
    opticalIsomers = 1,
    molecularWeight = (102.089,'amu'),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[1.01674,0.0773154,-0.000131084,1.28865e-07,-4.86782e-11,47938.8,33.0036], Tmin=(100,'K'), Tmax=(837.174,'K')), NASAPolynomial(coeffs=[2.9024,0.0421633,-2.12598e-05,4.10851e-09,-2.83957e-13,48539.2,29.713], Tmin=(837.174,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(397.787,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(286.849,'J/(mol*K)'), comment="""Thermo group additivity estimation: group(O2s-OsCs) + group(O2s-CsH) + group(O2s-OsH) + group(Cs-CsCsHH) + group(Cs-CsCsHH) + group(Cs-CsOsOsH) + group(Cs-CsHHH) + radical(CCOJ) + radical(RCCJ) + radical(RCCJC) + radical(Cs_P)"""),
)

species(
    label = '[CH2]CC[C]([O])O[O](1364)',
    structure = SMILES('[CH2]CC[C]([O])O[O]'),
    E0 = (355.345,'kJ/mol'),
    modes = [
        HarmonicOscillator(frequencies=([492.5,1135,1000,360,370,350,3000,3100,440,815,1455,1000,2750,2783.33,2816.67,2850,1425,1450,1225,1275,1270,1340,700,800,300,400,180,1652.36,1652.59],'cm^-1')),
        HinderedRotor(inertia=(0.189305,'amu*angstrom^2'), symmetry=1, barrier=(4.35249,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.189268,'amu*angstrom^2'), symmetry=1, barrier=(4.35163,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.18919,'amu*angstrom^2'), symmetry=1, barrier=(4.34984,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.189295,'amu*angstrom^2'), symmetry=1, barrier=(4.35226,'kJ/mol'), semiclassical=False),
    ],
    spinMultiplicity = 5,
    opticalIsomers = 1,
    molecularWeight = (102.089,'amu'),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[1.0684,0.0736602,-0.000116745,1.09917e-07,-4.0695e-11,42834.8,30.8425], Tmin=(100,'K'), Tmax=(826.424,'K')), NASAPolynomial(coeffs=[4.44447,0.0385788,-1.90555e-05,3.67186e-09,-2.54181e-13,42916.8,19.0705], Tmin=(826.424,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(355.345,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(291.007,'J/(mol*K)'), comment="""Thermo group additivity estimation: group(O2s-OsCs) + group(O2s-CsH) + group(O2s-OsH) + group(Cs-CsCsHH) + group(Cs-CsCsHH) + group(Cs-CsOsOsH) + group(Cs-CsHHH) + radical(Cs_P) + radical(RCCJ) + radical(CCOJ) + radical(ROOJ)"""),
)

species(
    label = '[CH2][CH][CH]C(O)O[O](8521)',
    structure = SMILES('[CH2][CH][CH]C(O)O[O]'),
    E0 = (319.255,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
    molecularWeight = (102.089,'amu'),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[0.916784,0.0776077,-0.000128196,1.18568e-07,-4.2209e-11,38499.1,34.8796], Tmin=(100,'K'), Tmax=(859.461,'K')), NASAPolynomial(coeffs=[5.45714,0.0350974,-1.66906e-05,3.13348e-09,-2.12455e-13,38508.3,18.2564], Tmin=(859.461,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(319.255,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(286.849,'J/(mol*K)'), comment="""Thermo group additivity estimation: group(O2s-OsCs) + group(O2s-CsH) + group(O2s-OsH) + group(Cs-CsCsHH) + group(Cs-CsCsHH) + group(Cs-CsOsOsH) + group(Cs-CsHHH) + radical(CCJCOOH) + radical(RCCJ) + radical(RCCJC) + radical(ROOJ)"""),
)

species(
    label = 'C[CH][CH]C([O])O[O](8279)',
    structure = SMILES('C[CH][CH]C([O])O[O]'),
    E0 = (339.714,'kJ/mol'),
    modes = [
        HarmonicOscillator(frequencies=([3000,3050,390,425,1340,1360,335,370,492.5,1135,1000,1380,1390,370,380,2900,435,2750,2800,2850,1350,1500,750,1050,1375,1000,180,2042.95,2043.05],'cm^-1')),
        HinderedRotor(inertia=(0.119499,'amu*angstrom^2'), symmetry=1, barrier=(2.74752,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.119333,'amu*angstrom^2'), symmetry=1, barrier=(2.74369,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.11943,'amu*angstrom^2'), symmetry=1, barrier=(2.74593,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.119418,'amu*angstrom^2'), symmetry=1, barrier=(2.74565,'kJ/mol'), semiclassical=False),
    ],
    spinMultiplicity = 5,
    opticalIsomers = 1,
    molecularWeight = (102.089,'amu'),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[1.22759,0.0743169,-0.000130828,1.32668e-07,-5.05575e-11,40945.1,32.8368], Tmin=(100,'K'), Tmax=(855.124,'K')), NASAPolynomial(coeffs=[0.953831,0.0438405,-2.16625e-05,4.13196e-09,-2.82578e-13,42153,40.9037], Tmin=(855.124,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(339.714,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(291.007,'J/(mol*K)'), comment="""Thermo group additivity estimation: group(O2s-OsCs) + group(O2s-CsH) + group(O2s-OsH) + group(Cs-CsCsHH) + group(Cs-CsCsHH) + group(Cs-CsOsOsH) + group(Cs-CsHHH) + radical(ROOJ) + radical(CCOJ) + radical(RCCJC) + radical(CCJCOOH)"""),
)

species(
    label = 'C[CH]C[C]([O])O[O](8280)',
    structure = SMILES('C[CH]C[C]([O])O[O]'),
    E0 = (344.545,'kJ/mol'),
    modes = [
        HarmonicOscillator(frequencies=([3025,407.5,1350,352.5,492.5,1135,1000,360,370,350,2750,2800,2850,1350,1500,750,1050,1375,1000,2750,2850,1437.5,1250,1305,750,350,180,2081.38,2082.88],'cm^-1')),
        HinderedRotor(inertia=(0.142297,'amu*angstrom^2'), symmetry=1, barrier=(3.27168,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.142267,'amu*angstrom^2'), symmetry=1, barrier=(3.27101,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.142352,'amu*angstrom^2'), symmetry=1, barrier=(3.27294,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.141767,'amu*angstrom^2'), symmetry=1, barrier=(3.25949,'kJ/mol'), semiclassical=False),
    ],
    spinMultiplicity = 5,
    opticalIsomers = 1,
    molecularWeight = (102.089,'amu'),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[1.24374,0.072974,-0.00012519,1.25511e-07,-4.76779e-11,41526.5,31.1619], Tmin=(100,'K'), Tmax=(851.686,'K')), NASAPolynomial(coeffs=[1.47381,0.042879,-2.10857e-05,4.02209e-09,-2.75423e-13,42539.6,36.2667], Tmin=(851.686,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(344.545,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(291.007,'J/(mol*K)'), comment="""Thermo group additivity estimation: group(O2s-OsCs) + group(O2s-CsH) + group(O2s-OsH) + group(Cs-CsCsHH) + group(Cs-CsCsHH) + group(Cs-CsOsOsH) + group(Cs-CsHHH) + radical(ROOJ) + radical(Cs_P) + radical(CCOJ) + radical(RCCJC)"""),
)

species(
    label = '[CH2][CH][CH]C([O])OO(8522)',
    structure = SMILES('[CH2][CH][CH]C([O])OO'),
    E0 = (392.955,'kJ/mol'),
    modes = [
        HarmonicOscillator(frequencies=([3000,3100,440,815,1455,1000,3615,1310,387.5,850,1000,1380,1390,370,380,2900,435,3000,3050,390,425,1340,1360,335,370,292.61,977.373,2340.88],'cm^-1')),
        HinderedRotor(inertia=(0.0659435,'amu*angstrom^2'), symmetry=1, barrier=(3.24539,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.0659435,'amu*angstrom^2'), symmetry=1, barrier=(3.24539,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.0659435,'amu*angstrom^2'), symmetry=1, barrier=(3.24539,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.0659435,'amu*angstrom^2'), symmetry=1, barrier=(3.24539,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.0659435,'amu*angstrom^2'), symmetry=1, barrier=(3.24539,'kJ/mol'), semiclassical=False),
    ],
    spinMultiplicity = 5,
    opticalIsomers = 1,
    molecularWeight = (102.089,'amu'),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[0.99958,0.0786698,-0.000136759,1.36059e-07,-5.1565e-11,47357.4,34.6822], Tmin=(100,'K'), Tmax=(841.969,'K')), NASAPolynomial(coeffs=[2.38987,0.0431116,-2.18289e-05,4.21653e-09,-2.90956e-13,48149.5,34.3084], Tmin=(841.969,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(392.955,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(286.849,'J/(mol*K)'), comment="""Thermo group additivity estimation: group(O2s-OsCs) + group(O2s-CsH) + group(O2s-OsH) + group(Cs-CsCsHH) + group(Cs-CsCsHH) + group(Cs-CsOsOsH) + group(Cs-CsHHH) + radical(CCJCOOH) + radical(RCCJC) + radical(CCOJ) + radical(RCCJ)"""),
)

species(
    label = '[CH2][CH]C[CH][O](749)',
    structure = SMILES('[CH2][CH]C[CH][O]'),
    E0 = (501.564,'kJ/mol'),
    modes = [
        HarmonicOscillator(frequencies=([2750,2850,1437.5,1250,1305,750,350,3000,3100,440,815,1455,1000,3000,3050,390,425,1340,1360,335,370,1976.59,1976.77,1977.46],'cm^-1')),
        HinderedRotor(inertia=(0.00270787,'amu*angstrom^2'), symmetry=1, barrier=(7.50732,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.326656,'amu*angstrom^2'), symmetry=1, barrier=(7.51047,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(2.47397,'amu*angstrom^2'), symmetry=1, barrier=(56.8814,'kJ/mol'), semiclassical=False),
    ],
    spinMultiplicity = 5,
    opticalIsomers = 1,
    molecularWeight = (70.0898,'amu'),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[2.0673,0.0533416,-9.28541e-05,9.68591e-08,-3.74643e-11,60383.4,24.175], Tmin=(100,'K'), Tmax=(867.002,'K')), NASAPolynomial(coeffs=[-0.0667086,0.0364534,-1.73841e-05,3.26323e-09,-2.20956e-13,61758.2,39.9603], Tmin=(867.002,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(501.564,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(245.277,'J/(mol*K)'), comment="""Thermo group additivity estimation: group(O2s-CsH) + group(Cs-CsCsHH) + group(Cs-CsCsHH) + group(Cs-CsOsHH) + group(Cs-CsHHH) + radical(RCCJ) + radical(RCCJC) + radical(CCOJ) + radical(CCsJOH)"""),
)

species(
    label = '[CH2][CH][CH]C([O])O[O](8523)',
    structure = SMILES('[CH2][CH][CH]C([O])O[O]'),
    E0 = (544.96,'kJ/mol'),
    modes = [
        HarmonicOscillator(frequencies=([3000,3100,440,815,1455,1000,492.5,1135,1000,1380,1390,370,380,2900,435,3000,3050,390,425,1340,1360,335,370,180,1763.33,1767.92],'cm^-1')),
        HinderedRotor(inertia=(0.137553,'amu*angstrom^2'), symmetry=1, barrier=(3.16262,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.136369,'amu*angstrom^2'), symmetry=1, barrier=(3.13539,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.0014424,'amu*angstrom^2'), symmetry=1, barrier=(3.18117,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.139031,'amu*angstrom^2'), symmetry=1, barrier=(3.19659,'kJ/mol'), semiclassical=False),
    ],
    spinMultiplicity = 6,
    opticalIsomers = 1,
    molecularWeight = (101.081,'amu'),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[1.25748,0.0751339,-0.000141545,1.45194e-07,-5.5062e-11,65628.1,34.5161], Tmin=(100,'K'), Tmax=(863.471,'K')), NASAPolynomial(coeffs=[0.797161,0.0414594,-2.0843e-05,3.97709e-09,-2.7094e-13,67042.4,44.3989], Tmin=(863.471,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(544.96,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(266.063,'J/(mol*K)'), comment="""Thermo group additivity estimation: group(O2s-OsCs) + group(O2s-CsH) + group(O2s-OsH) + group(Cs-CsCsHH) + group(Cs-CsCsHH) + group(Cs-CsOsOsH) + group(Cs-CsHHH) + radical(RCCJC) + radical(RCCJ) + radical(CCJCOOH) + radical(ROOJ) + radical(CCOJ)"""),
)

species(
    label = '[CH2][CH]C[C]([O])O[O](8524)',
    structure = SMILES('[CH2][CH]C[C]([O])O[O]'),
    E0 = (549.791,'kJ/mol'),
    modes = [
        HarmonicOscillator(frequencies=([3025,407.5,1350,352.5,492.5,1135,1000,360,370,350,3000,3100,440,815,1455,1000,2750,2850,1437.5,1250,1305,750,350,1786.05,1786.4,1786.52],'cm^-1')),
        HinderedRotor(inertia=(0.159516,'amu*angstrom^2'), symmetry=1, barrier=(3.66759,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.15941,'amu*angstrom^2'), symmetry=1, barrier=(3.66514,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.159384,'amu*angstrom^2'), symmetry=1, barrier=(3.66455,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.159307,'amu*angstrom^2'), symmetry=1, barrier=(3.66277,'kJ/mol'), semiclassical=False),
    ],
    spinMultiplicity = 6,
    opticalIsomers = 1,
    molecularWeight = (101.081,'amu'),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[1.27263,0.0738032,-0.000135949,1.38091e-07,-5.22038e-11,66209.5,32.8448], Tmin=(100,'K'), Tmax=(860.907,'K')), NASAPolynomial(coeffs=[1.32211,0.0404893,-2.02611e-05,3.86598e-09,-2.63681e-13,67427.1,39.7341], Tmin=(860.907,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(549.791,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(266.063,'J/(mol*K)'), comment="""Thermo group additivity estimation: group(O2s-OsCs) + group(O2s-CsH) + group(O2s-OsH) + group(Cs-CsCsHH) + group(Cs-CsCsHH) + group(Cs-CsOsOsH) + group(Cs-CsHHH) + radical(Cs_P) + radical(ROOJ) + radical(RCCJC) + radical(CCOJ) + radical(RCCJ)"""),
)

species(
    label = 'C=C[CH]C(O)O[O](8525)',
    structure = SMILES('C=C[CH]C(O)O[O]'),
    E0 = (-36.1539,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
    molecularWeight = (102.089,'amu'),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[0.955637,0.0682028,-7.07217e-05,3.87456e-08,-8.62659e-12,-4239.67,25.5939], Tmin=(100,'K'), Tmax=(1077.22,'K')), NASAPolynomial(coeffs=[12.3556,0.0258709,-1.17746e-05,2.2639e-09,-1.59794e-13,-6695.68,-30.2522], Tmin=(1077.22,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(-36.1539,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(291.007,'J/(mol*K)'), comment="""Thermo group additivity estimation: group(O2s-OsCs) + group(O2s-CsH) + group(O2s-OsH) + group(Cs-(Cds-Cds)CsHH) + group(Cs-CsOsOsH) + group(Cds-CdsCsH) + group(Cds-CdsHH) + radical(C=CCJCO) + radical(ROOJ)"""),
)

species(
    label = 'C=CC[C]([O])OO(1285)',
    structure = SMILES('C=CC[C]([O])OO'),
    E0 = (125.876,'kJ/mol'),
    modes = [
        HarmonicOscillator(frequencies=([3615,1310,387.5,850,1000,360,370,350,2950,3100,1380,975,1025,1650,3010,987.5,1337.5,450,1655,2750,2850,1437.5,1250,1305,750,350,259.605,259.606,1846.41],'cm^-1')),
        HinderedRotor(inertia=(0.213021,'amu*angstrom^2'), symmetry=1, barrier=(10.1877,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.00421111,'amu*angstrom^2'), symmetry=1, barrier=(10.1877,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.695368,'amu*angstrom^2'), symmetry=1, barrier=(33.2558,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.695363,'amu*angstrom^2'), symmetry=1, barrier=(33.2558,'kJ/mol'), semiclassical=False),
    ],
    spinMultiplicity = 3,
    opticalIsomers = 1,
    molecularWeight = (102.089,'amu'),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[1.58912,0.0611461,-5.58922e-05,1.10379e-08,1.51329e-11,15218.3,26.9869], Tmin=(100,'K'), Tmax=(527.595,'K')), NASAPolynomial(coeffs=[6.24814,0.0369151,-1.8536e-05,3.68192e-09,-2.63025e-13,14572.3,6.02591], Tmin=(527.595,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(125.876,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(291.007,'J/(mol*K)'), comment="""Thermo group additivity estimation: group(O2s-OsCs) + group(O2s-CsH) + group(O2s-OsH) + group(Cs-(Cds-Cds)CsHH) + group(Cs-CsOsOsH) + group(Cds-CdsCsH) + group(Cds-CdsHH) + radical(Cs_P) + radical(CCOJ)"""),
)

species(
    label = '[CH2]CCC(=O)O[O](8526)',
    structure = SMILES('[CH2]CCC(=O)O[O]'),
    E0 = (-95.3158,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
    molecularWeight = (102.089,'amu'),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[1.07985,0.0566216,-4.36056e-05,1.71334e-08,-2.71077e-12,-11352.4,28.2118], Tmin=(100,'K'), Tmax=(1495.64,'K')), NASAPolynomial(coeffs=[14.132,0.0217143,-8.59662e-06,1.52848e-09,-1.02379e-13,-15256.7,-40.0117], Tmin=(1495.64,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(-95.3158,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(291.007,'J/(mol*K)'), comment="""Thermo group additivity estimation: group(O2s-O2s(Cds-O2d)) + group(O2s-OsH) + group(Cs-CsCsHH) + group(Cs-(Cds-O2d)CsHH) + group(Cs-CsHHH) + group(Cds-OdCsOs) + radical(RCCJ) + radical(C(=O)OOJ)"""),
)

species(
    label = 'C=C[CH]C([O])OO(8527)',
    structure = SMILES('C=C[CH]C([O])OO'),
    E0 = (37.5465,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
    molecularWeight = (102.089,'amu'),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[1.26094,0.0663691,-6.78518e-05,3.94626e-08,-9.91046e-12,4609.3,24.6161], Tmin=(100,'K'), Tmax=(929.838,'K')), NASAPolynomial(coeffs=[8.47819,0.0353218,-1.77667e-05,3.5531e-09,-2.55677e-13,3267.12,-9.67798], Tmin=(929.838,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(37.5465,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(291.007,'J/(mol*K)'), comment="""Thermo group additivity estimation: group(O2s-OsCs) + group(O2s-CsH) + group(O2s-OsH) + group(Cs-(Cds-Cds)CsHH) + group(Cs-CsOsOsH) + group(Cds-CdsCsH) + group(Cds-CdsHH) + radical(CCOJ) + radical(C=CCJCO)"""),
)

species(
    label = 'O2(S)(5486)',
    structure = SMILES('O=O'),
    E0 = (85.6848,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
    molecularWeight = (31.9988,'amu'),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[3.53732,-0.00121571,5.31618e-06,-4.89443e-09,1.45845e-12,10304.5,4.68368], Tmin=(100,'K'), Tmax=(1074.56,'K')), NASAPolynomial(coeffs=[3.15382,0.00167804,-7.69971e-07,1.51275e-10,-1.08782e-14,10302.3,6.16754], Tmin=(1074.56,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(85.6848,'kJ/mol'), Cp0=(29.1007,'J/(mol*K)'), CpInf=(37.4151,'J/(mol*K)'), label="""O2(S)""", comment="""Thermo library: primaryThermoLibrary"""),
)

species(
    label = 'O(T)(63)',
    structure = SMILES('[O]'),
    E0 = (243.034,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
    molecularWeight = (15.9994,'amu'),
    collisionModel = TransportData(shapeIndex=0, epsilon=(665.16,'J/mol'), sigma=(2.75,'angstroms'), dipoleMoment=(0,'C*m'), polarizability=(0,'angstroms^3'), rotrelaxcollnum=0.0, comment="""GRI-Mech"""),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[2.5,9.24385e-15,-1.3678e-17,6.66185e-21,-1.00107e-24,29230.2,4.09104], Tmin=(100,'K'), Tmax=(3459.6,'K')), NASAPolynomial(coeffs=[2.5,9.20456e-12,-3.58608e-15,6.15199e-19,-3.92042e-23,29230.2,4.09104], Tmin=(3459.6,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(243.034,'kJ/mol'), Cp0=(20.7862,'J/(mol*K)'), CpInf=(20.7862,'J/(mol*K)'), label="""O(T)""", comment="""Thermo library: primaryThermoLibrary"""),
)

species(
    label = '[CH2][CH]CC1OO1(8528)',
    structure = SMILES('[CH2][CH]CC1OO1'),
    E0 = (289.541,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
    molecularWeight = (86.0892,'amu'),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[1.9261,0.0373484,-1.65452e-06,-2.6514e-08,1.39465e-11,34906,24.9814], Tmin=(100,'K'), Tmax=(927.472,'K')), NASAPolynomial(coeffs=[11.1308,0.0187514,-5.7043e-06,9.27202e-10,-6.23588e-14,32291,-23.6256], Tmin=(927.472,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(289.541,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(270.22,'J/(mol*K)'), comment="""Thermo group additivity estimation: group(O2s-OsCs) + group(O2s-OsCs) + group(Cs-CsCsHH) + group(Cs-CsCsHH) + group(Cs-CsOsOsH) + group(Cs-CsHHH) + ring(dioxirane) + radical(RCCJC) + radical(RCCJ)"""),
)

species(
    label = '[CH2]C1CC([O])O1(712)',
    structure = SMILES('[CH2]C1CC([O])O1'),
    E0 = (103.242,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
    molecularWeight = (86.0892,'amu'),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[1.52144,0.0463373,-3.51214e-05,1.57586e-08,-2.87925e-12,12513.3,21.3473], Tmin=(100,'K'), Tmax=(1468.21,'K')), NASAPolynomial(coeffs=[9.29478,0.0214877,-5.98235e-06,8.24147e-10,-4.62333e-14,10626.5,-17.7921], Tmin=(1468.21,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(103.242,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(278.535,'J/(mol*K)'), comment="""Thermo group additivity estimation: group(O2s-CsCs) + group(O2s-CsH) + group(Cs-CsCsOsH) + group(Cs-CsCsHH) + group(Cs-CsOsOsH) + group(Cs-CsHHH) + ring(Oxetane) + radical(CCOJ) + radical(CJC(C)OC)"""),
)

species(
    label = '[CH2]C([CH2])C([O])O[O](8175)',
    structure = SMILES('[CH2]C([CH2])C([O])O[O]'),
    E0 = (345.975,'kJ/mol'),
    modes = [
        HarmonicOscillator(frequencies=([492.5,1135,1000,1380,1383.33,1386.67,1390,370,373.333,376.667,380,2800,3000,430,440,3000,3033.33,3066.67,3100,415,465,780,850,1435,1475,900,1100,1279.19,1279.35],'cm^-1')),
        HinderedRotor(inertia=(0.175661,'amu*angstrom^2'), symmetry=1, barrier=(4.03879,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.175724,'amu*angstrom^2'), symmetry=1, barrier=(4.04025,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.175436,'amu*angstrom^2'), symmetry=1, barrier=(4.03361,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.175448,'amu*angstrom^2'), symmetry=1, barrier=(4.03389,'kJ/mol'), semiclassical=False),
    ],
    spinMultiplicity = 5,
    opticalIsomers = 1,
    molecularWeight = (102.089,'amu'),
    collisionModel = TransportData(shapeIndex=2, epsilon=(4236.76,'J/mol'), sigma=(7.09873,'angstroms'), dipoleMoment=(0,'C*m'), polarizability=(0,'angstroms^3'), rotrelaxcollnum=0, comment="""Epsilon & sigma estimated with Tc=661.77 K, Pc=26.87 bar (from Joback method)"""),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[1.01052,0.0733357,-0.000111213,9.73465e-08,-3.31066e-11,41711.7,31.1696], Tmin=(100,'K'), Tmax=(885.58,'K')), NASAPolynomial(coeffs=[6.04596,0.0334035,-1.44622e-05,2.59603e-09,-1.71385e-13,41493.8,11.2937], Tmin=(885.58,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(345.975,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(291.007,'J/(mol*K)'), comment="""Thermo group additivity estimation: group(O2s-OsCs) + group(O2s-CsH) + group(O2s-OsH) + group(Cs-CsCsCsH) + group(Cs-CsOsOsH) + group(Cs-CsHHH) + group(Cs-CsHHH) + radical(Isobutyl) + radical(CCOJ) + radical(Isobutyl) + radical(ROOJ)"""),
)

species(
    label = '[O]O(16)',
    structure = SMILES('[O]O'),
    E0 = (-8.19602,'kJ/mol'),
    modes = [
        HarmonicOscillator(frequencies=([1036.72,2034.11,2034.11],'cm^-1')),
    ],
    spinMultiplicity = 2,
    opticalIsomers = 1,
    molecularWeight = (33.0067,'amu'),
    collisionModel = TransportData(shapeIndex=2, epsilon=(892.977,'J/mol'), sigma=(3.458,'angstroms'), dipoleMoment=(0,'C*m'), polarizability=(0,'angstroms^3'), rotrelaxcollnum=1.0, comment="""GRI-Mech"""),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[4.04595,-0.00173474,1.0377e-05,-1.02207e-08,3.3493e-12,-986.755,4.63579], Tmin=(100,'K'), Tmax=(932.129,'K')), NASAPolynomial(coeffs=[3.21022,0.00367946,-1.27704e-06,2.18051e-10,-1.46343e-14,-910.359,8.18305], Tmin=(932.129,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(-8.19602,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(58.2013,'J/(mol*K)'), comment="""Thermo group additivity estimation: group(O2s-OsH) + group(O2s-OsH) + radical(HOOJ)"""),
)

species(
    label = '[CH2][CH]C=C[O](5307)',
    structure = SMILES('[CH2][CH]C=C[O]'),
    E0 = (262.507,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
    molecularWeight = (69.0819,'amu'),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[2.26611,0.0245108,2.84591e-05,-6.02061e-08,2.6419e-11,31647.2,18.9439], Tmin=(100,'K'), Tmax=(940.567,'K')), NASAPolynomial(coeffs=[13.977,0.00916547,-2.0216e-06,3.48798e-10,-2.92061e-14,27920.1,-44.9392], Tmin=(940.567,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(262.507,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(224.491,'J/(mol*K)'), comment="""Thermo group additivity estimation: group(O2s-(Cds-Cd)H) + group(Cs-(Cds-Cds)CsHH) + group(Cs-CsHHH) + group(Cds-CdsCsH) + group(Cds-CdsOsH) + radical(RCCJ) + radical(C=COJ) + radical(Allyl_S)"""),
)

species(
    label = '[CH2][CH]CC1OOO1(8529)',
    structure = SMILES('[CH2][CH]CC1OOO1'),
    E0 = (330.754,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
    molecularWeight = (102.089,'amu'),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[1.82278,0.0334429,2.68064e-05,-5.62529e-08,2.27059e-11,39871.5,28.4501], Tmin=(100,'K'), Tmax=(1009.05,'K')), NASAPolynomial(coeffs=[13.1292,0.023301,-9.66765e-06,1.90348e-09,-1.40912e-13,35824.3,-34.947], Tmin=(1009.05,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(330.754,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(295.164,'J/(mol*K)'), comment="""Thermo group additivity estimation: group(O2s-OsCs) + group(O2s-OsCs) + group(O2s-OsOs) + group(Cs-CsCsHH) + group(Cs-CsCsHH) + group(Cs-CsOsOsH) + group(Cs-CsHHH) + ring(Cyclobutane) + radical(RCCJC) + radical(RCCJ)"""),
)

species(
    label = '[CH2]C1CC(O[O])O1(8177)',
    structure = SMILES('[CH2]C1CC(O[O])O1'),
    E0 = (101.046,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
    molecularWeight = (102.089,'amu'),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[0.710022,0.0629746,-6.11704e-05,3.24116e-08,-6.64336e-12,12279.5,24.9125], Tmin=(100,'K'), Tmax=(1356.08,'K')), NASAPolynomial(coeffs=[13.1804,0.0188715,-4.29036e-06,4.6841e-10,-2.07063e-14,9570.38,-36.5669], Tmin=(1356.08,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(101.046,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(299.321,'J/(mol*K)'), comment="""Thermo group additivity estimation: group(O2s-CsCs) + group(O2s-OsCs) + group(O2s-OsH) + group(Cs-CsCsOsH) + group(Cs-CsCsHH) + group(Cs-CsOsOsH) + group(Cs-CsHHH) + ring(Oxetane) + radical(ROOJ) + radical(CJC(C)OC)"""),
)

species(
    label = '[CH2]C1CC([O])OO1(8530)',
    structure = SMILES('[CH2]C1CC([O])OO1'),
    E0 = (78.5931,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
    molecularWeight = (102.089,'amu'),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[1.64503,0.041012,7.8691e-06,-4.43437e-08,2.22958e-11,9547.45,23.9347], Tmin=(100,'K'), Tmax=(902.404,'K')), NASAPolynomial(coeffs=[12.9918,0.019707,-4.90655e-06,6.95266e-10,-4.44424e-14,6319.16,-36.1822], Tmin=(902.404,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(78.5931,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(303.478,'J/(mol*K)'), comment="""Thermo group additivity estimation: group(O2s-OsCs) + group(O2s-OsCs) + group(O2s-CsH) + group(Cs-CsCsOsH) + group(Cs-CsCsHH) + group(Cs-CsOsOsH) + group(Cs-CsHHH) + ring(12dioxolane) + radical(CCOJ) + radical(CJCOOH)"""),
)

species(
    label = 'C=CCC([O])O[O](8172)',
    structure = SMILES('C=CCC([O])O[O]'),
    E0 = (72.6346,'kJ/mol'),
    modes = [
        HarmonicOscillator(frequencies=([492.5,1135,1000,2950,3100,1380,975,1025,1650,3010,987.5,1337.5,450,1655,1380,1390,370,380,2900,435,2750,2850,1437.5,1250,1305,750,350,182.993,185.347,3616.93],'cm^-1')),
        HinderedRotor(inertia=(0.231258,'amu*angstrom^2'), symmetry=1, barrier=(5.59521,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(1.4956,'amu*angstrom^2'), symmetry=1, barrier=(36.8098,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.236773,'amu*angstrom^2'), symmetry=1, barrier=(5.59053,'kJ/mol'), semiclassical=False),
    ],
    spinMultiplicity = 3,
    opticalIsomers = 1,
    molecularWeight = (102.089,'amu'),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[1.3619,0.0636604,-8.35206e-05,7.22424e-08,-2.65514e-11,8825.5,27.7916], Tmin=(100,'K'), Tmax=(759.251,'K')), NASAPolynomial(coeffs=[5.23389,0.0368469,-1.78746e-05,3.474e-09,-2.44003e-13,8422.43,11.3955], Tmin=(759.251,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(72.6346,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(295.164,'J/(mol*K)'), comment="""Thermo group additivity estimation: group(O2s-OsCs) + group(O2s-CsH) + group(O2s-OsH) + group(Cs-(Cds-Cds)CsHH) + group(Cs-CsOsOsH) + group(Cds-CdsCsH) + group(Cds-CdsHH) + radical(ROOJ) + radical(CCOJ)"""),
)

species(
    label = '[CH2][CH]CC([O])[O](2311)',
    structure = SMILES('[CH2][CH]CC([O])[O]'),
    E0 = (346.74,'kJ/mol'),
    modes = [
        HarmonicOscillator(frequencies=([3000,3100,440,815,1455,1000,2750,2850,1437.5,1250,1305,750,350,1380,1390,370,380,2900,435,3025,407.5,1350,352.5,180,1796.73,1796.74,1796.82],'cm^-1')),
        HinderedRotor(inertia=(0.110343,'amu*angstrom^2'), symmetry=1, barrier=(2.53699,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.110447,'amu*angstrom^2'), symmetry=1, barrier=(2.5394,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.00110859,'amu*angstrom^2'), symmetry=1, barrier=(2.53965,'kJ/mol'), semiclassical=False),
    ],
    spinMultiplicity = 5,
    opticalIsomers = 1,
    molecularWeight = (86.0892,'amu'),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[2.01876,0.056897,-0.000101745,1.13239e-07,-4.62652e-11,41761.7,28.1242], Tmin=(100,'K'), Tmax=(846.772,'K')), NASAPolynomial(coeffs=[-2.81621,0.0461229,-2.31154e-05,4.45346e-09,-3.06969e-13,43785.6,57.7617], Tmin=(846.772,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(346.74,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(270.22,'J/(mol*K)'), comment="""Thermo group additivity estimation: group(O2s-CsH) + group(O2s-CsH) + group(Cs-CsCsHH) + group(Cs-CsCsHH) + group(Cs-CsOsOsH) + group(Cs-CsHHH) + radical(CCOJ) + radical(RCCJ) + radical(CCOJ) + radical(RCCJC)"""),
)

species(
    label = '[CH2]C([O])O[O](607)',
    structure = SMILES('[CH2]C([O])O[O]'),
    E0 = (206.375,'kJ/mol'),
    modes = [
        HarmonicOscillator(frequencies=([492.5,1135,1000,1380,1390,370,380,2900,435,3000,3100,440,815,1455,1000,1967.16],'cm^-1')),
        HinderedRotor(inertia=(0.24983,'amu*angstrom^2'), symmetry=1, barrier=(5.74408,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.249705,'amu*angstrom^2'), symmetry=1, barrier=(5.74121,'kJ/mol'), semiclassical=False),
    ],
    spinMultiplicity = 4,
    opticalIsomers = 1,
    molecularWeight = (75.0434,'amu'),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[2.29491,0.0444229,-7.76836e-05,7.63138e-08,-2.8547e-11,24875.9,21.4728], Tmin=(100,'K'), Tmax=(845.116,'K')), NASAPolynomial(coeffs=[3.56685,0.022769,-1.15016e-05,2.21713e-09,-1.52683e-13,25219.2,18.8536], Tmin=(845.116,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(206.375,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(174.604,'J/(mol*K)'), comment="""Thermo group additivity estimation: group(O2s-OsCs) + group(O2s-CsH) + group(O2s-OsH) + group(Cs-CsOsOsH) + group(Cs-CsHHH) + radical(ROOJ) + radical(CCOJ) + radical(CJCOOH)"""),
)

species(
    label = '[CH][CH2](721)',
    structure = SMILES('[CH][CH2]'),
    E0 = (556.742,'kJ/mol'),
    modes = [
        HarmonicOscillator(frequencies=([3000,3100,440,815,1455,1000,1101.59,1101.66],'cm^-1')),
        HinderedRotor(inertia=(0.00420677,'amu*angstrom^2'), symmetry=1, barrier=(3.62356,'kJ/mol'), semiclassical=False),
    ],
    spinMultiplicity = 4,
    opticalIsomers = 1,
    molecularWeight = (27.0452,'amu'),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[3.77493,-0.000462567,3.18167e-05,-4.30783e-08,1.77606e-11,66973.8,8.79001], Tmin=(100,'K'), Tmax=(870.354,'K')), NASAPolynomial(coeffs=[6.06996,0.00332438,5.85464e-07,-2.32999e-10,1.82455e-14,66031.4,-5.08252], Tmin=(870.354,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(556.742,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(103.931,'J/(mol*K)'), comment="""Thermo group additivity estimation: group(Cs-CsHHH) + group(Cs-CsHHH) + radical(CCJ2_triplet) + radical(CCJ)"""),
)

species(
    label = '[CH2][CH]C[CH]O[O](8531)',
    structure = SMILES('[CH2][CH]C[CH]O[O]'),
    E0 = (507.653,'kJ/mol'),
    modes = [
        HarmonicOscillator(frequencies=([3000,3050,390,425,1340,1360,335,370,492.5,1135,1000,3000,3100,440,815,1455,1000,2750,2850,1437.5,1250,1305,750,350,1436.17,1437.8],'cm^-1')),
        HinderedRotor(inertia=(0.161881,'amu*angstrom^2'), symmetry=1, barrier=(3.72197,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.160856,'amu*angstrom^2'), symmetry=1, barrier=(3.69839,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.160181,'amu*angstrom^2'), symmetry=1, barrier=(3.68287,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.160357,'amu*angstrom^2'), symmetry=1, barrier=(3.68693,'kJ/mol'), semiclassical=False),
    ],
    spinMultiplicity = 5,
    opticalIsomers = 1,
    molecularWeight = (86.0892,'amu'),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[1.52198,0.0639485,-0.000103896,9.96062e-08,-3.65091e-11,61136.8,27.3928], Tmin=(100,'K'), Tmax=(866.778,'K')), NASAPolynomial(coeffs=[2.91381,0.0354909,-1.65172e-05,3.07179e-09,-2.07174e-13,61723.2,25.6518], Tmin=(866.778,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(507.653,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(266.063,'J/(mol*K)'), comment="""Thermo group additivity estimation: group(O2s-OsCs) + group(O2s-OsH) + group(Cs-CsCsHH) + group(Cs-CsCsHH) + group(Cs-CsOsHH) + group(Cs-CsHHH) + radical(RCCJ) + radical(ROOJ) + radical(CCsJOOH) + radical(RCCJC)"""),
)

species(
    label = '[CH2][C]CC([O])O[O](8532)',
    structure = SMILES('[CH2][C]CC([O])O[O]'),
    E0 = (598.314,'kJ/mol'),
    modes = [
        HarmonicOscillator(frequencies=([3000,3100,440,815,1455,1000,492.5,1135,1000,1380,1390,370,380,2900,435,2750,2850,1437.5,1250,1305,750,350,216.276,223.6,1654.16,1655.18],'cm^-1')),
        HinderedRotor(inertia=(0.172245,'amu*angstrom^2'), symmetry=1, barrier=(5.81104,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.173677,'amu*angstrom^2'), symmetry=1, barrier=(5.81035,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.173346,'amu*angstrom^2'), symmetry=1, barrier=(5.80764,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(1.43128,'amu*angstrom^2'), symmetry=1, barrier=(48.7998,'kJ/mol'), semiclassical=False),
    ],
    spinMultiplicity = 6,
    opticalIsomers = 1,
    molecularWeight = (101.081,'amu'),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[0.966791,0.0761155,-0.000126217,1.17729e-07,-4.27061e-11,72060.7,30.4251], Tmin=(100,'K'), Tmax=(834.744,'K')), NASAPolynomial(coeffs=[5.62959,0.0346521,-1.73518e-05,3.34416e-09,-2.30844e-13,71948.3,12.7621], Tmin=(834.744,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(598.314,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(266.063,'J/(mol*K)'), comment="""Thermo group additivity estimation: group(O2s-OsCs) + group(O2s-CsH) + group(O2s-OsH) + group(Cs-CsCsHH) + group(Cs-CsCsHH) + group(Cs-CsOsOsH) + group(Cs-CsHHH) + radical(RCCJ) + radical(ROOJ) + radical(CCOJ) + radical(CCJ2_triplet)"""),
)

species(
    label = '[CH][CH]CC([O])O[O](8533)',
    structure = SMILES('[CH][CH]CC([O])O[O]'),
    E0 = (587.514,'kJ/mol'),
    modes = [
        HarmonicOscillator(frequencies=([3025,407.5,1350,352.5,492.5,1135,1000,1380,1390,370,380,2900,435,2750,2850,1437.5,1250,1305,750,350,208.909,901.275,1080.99,1297.19,1577.23,1876.45],'cm^-1')),
        HinderedRotor(inertia=(0.122036,'amu*angstrom^2'), symmetry=1, barrier=(3.27876,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.122036,'amu*angstrom^2'), symmetry=1, barrier=(3.27876,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.122036,'amu*angstrom^2'), symmetry=1, barrier=(3.27876,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.122036,'amu*angstrom^2'), symmetry=1, barrier=(3.27876,'kJ/mol'), semiclassical=False),
    ],
    spinMultiplicity = 6,
    opticalIsomers = 1,
    molecularWeight = (101.081,'amu'),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[1.14541,0.0753907,-0.000134531,1.33173e-07,-4.96422e-11,70752.2,31.8313], Tmin=(100,'K'), Tmax=(857.857,'K')), NASAPolynomial(coeffs=[2.63873,0.038988,-1.94031e-05,3.69946e-09,-2.52513e-13,71579.3,31.1696], Tmin=(857.857,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(587.514,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(266.063,'J/(mol*K)'), comment="""Thermo group additivity estimation: group(O2s-OsCs) + group(O2s-CsH) + group(O2s-OsH) + group(Cs-CsCsHH) + group(Cs-CsCsHH) + group(Cs-CsOsOsH) + group(Cs-CsHHH) + radical(CCJ2_triplet) + radical(RCCJC) + radical(ROOJ) + radical(CCOJ)"""),
)

species(
    label = 'N2',
    structure = SMILES('N#N'),
    E0 = (-8.64289,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
    molecularWeight = (28.0135,'amu'),
    collisionModel = TransportData(shapeIndex=1, epsilon=(810.913,'J/mol'), sigma=(3.621,'angstroms'), dipoleMoment=(0,'C*m'), polarizability=(1.76,'angstroms^3'), rotrelaxcollnum=4.0, comment="""GRI-Mech"""),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[3.53101,-0.000123661,-5.02999e-07,2.43531e-09,-1.40881e-12,-1046.98,2.96747], Tmin=(200,'K'), Tmax=(1000,'K')), NASAPolynomial(coeffs=[2.95258,0.0013969,-4.92632e-07,7.8601e-11,-4.60755e-15,-923.949,5.87189], Tmin=(1000,'K'), Tmax=(6000,'K'))], Tmin=(200,'K'), Tmax=(6000,'K'), E0=(-8.64289,'kJ/mol'), Cp0=(29.1007,'J/(mol*K)'), CpInf=(37.4151,'J/(mol*K)'), label="""N2""", comment="""Thermo library: primaryThermoLibrary"""),
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

transitionState(
    label = 'TS1',
    E0 = (344.545,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS2',
    E0 = (409.13,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS3',
    E0 = (536.191,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS4',
    E0 = (422.534,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS5',
    E0 = (345.952,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS6',
    E0 = (344.545,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS7',
    E0 = (487.517,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS8',
    E0 = (502.581,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS9',
    E0 = (532.959,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS10',
    E0 = (471.481,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS11',
    E0 = (419.826,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS12',
    E0 = (488.875,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS13',
    E0 = (473.013,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS14',
    E0 = (476.844,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS15',
    E0 = (492.937,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS16',
    E0 = (713.18,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS17',
    E0 = (756.765,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS18',
    E0 = (761.596,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS19',
    E0 = (407.945,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS20',
    E0 = (407.945,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS21',
    E0 = (407.945,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS22',
    E0 = (369.518,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS23',
    E0 = (344.545,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS24',
    E0 = (534.075,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS25',
    E0 = (427.443,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS26',
    E0 = (503.294,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS27',
    E0 = (483.914,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS28',
    E0 = (352.829,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS29',
    E0 = (352.829,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS30',
    E0 = (351.658,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS31',
    E0 = (344.545,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS32',
    E0 = (594.941,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS33',
    E0 = (797.433,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS34',
    E0 = (914.534,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS35',
    E0 = (810.118,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS36',
    E0 = (799.318,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

reaction(
    label = 'reaction1',
    reactants = ['[CH2][CH]CC([O])O[O](8173)'],
    products = ['[CH2]C=C(87)', '[O]OC=O(5472)'],
    transitionState = 'TS1',
    kinetics = Arrhenius(A=(5e+12,'s^-1'), n=0, Ea=(0,'kJ/mol'), T0=(1,'K'), Tmin=(300,'K'), Tmax=(1500,'K'), comment="""Exact match found for rate rule [RJJ]
Euclidian distance = 0
family: 1,4_Linear_birad_scission"""),
)

reaction(
    label = 'reaction2',
    reactants = ['H(8)', 'C=C[CH]C([O])O[O](8518)'],
    products = ['[CH2][CH]CC([O])O[O](8173)'],
    transitionState = 'TS2',
    kinetics = Arrhenius(A=(8.22e+08,'cm^3/(mol*s)'), n=1.533, Ea=(7.77387,'kJ/mol'), T0=(1,'K'), comment="""From training reaction 192 used for Cd_R;HJ
Exact match found for rate rule [Cd_R;HJ]
Euclidian distance = 0
family: R_Addition_MultipleBond"""),
)

reaction(
    label = 'reaction3',
    reactants = ['C=CC[C]([O])O[O](5392)', 'H(8)'],
    products = ['[CH2][CH]CC([O])O[O](8173)'],
    transitionState = 'TS3',
    kinetics = Arrhenius(A=(1.83701,'m^3/(mol*s)'), n=1.71338, Ea=(46.5052,'kJ/mol'), T0=(1,'K'), comment="""Estimated using average of templates [CO_O;HJ] + [CO-NdNd_O;YJ] for rate rule [CO-NdNd_O;HJ]
Euclidian distance = 1.0
family: R_Addition_MultipleBond"""),
)

reaction(
    label = 'reaction4',
    reactants = ['[CH2]C=C(87)', '[O][CH]O[O](8201)'],
    products = ['[CH2][CH]CC([O])O[O](8173)'],
    transitionState = 'TS4',
    kinetics = Arrhenius(A=(0.246938,'m^3/(mol*s)'), n=2.00579, Ea=(36.0234,'kJ/mol'), T0=(1,'K'), comment="""Estimated using an average for rate rule [Cd_R;CJ]
Euclidian distance = 0
Multiplied by reaction path degeneracy 2.0
family: R_Addition_MultipleBond"""),
)

reaction(
    label = 'reaction5',
    reactants = ['[O]OC=O(5472)', '[CH2][CH][CH2](6136)'],
    products = ['[CH2][CH]CC([O])O[O](8173)'],
    transitionState = 'TS5',
    kinetics = Arrhenius(A=(13.7458,'m^3/(mol*s)'), n=1.39198, Ea=(57.1551,'kJ/mol'), T0=(1,'K'), comment="""Estimated using average of templates [CO_O;CJ] + [CO-NdH_O;YJ] for rate rule [CO-NdH_O;CJ]
Euclidian distance = 1.0
Multiplied by reaction path degeneracy 2.0
family: R_Addition_MultipleBond"""),
)

reaction(
    label = 'reaction6',
    reactants = ['O2(2)', 'C=CC[CH][O](692)'],
    products = ['[CH2][CH]CC([O])O[O](8173)'],
    transitionState = 'TS6',
    kinetics = Arrhenius(A=(0.0849,'cm^3/(mol*s)'), n=3.486, Ea=(123.518,'kJ/mol'), T0=(1,'K'), Tmin=(300,'K'), Tmax=(1600,'K'), comment="""Estimated using template [CO-CsH_O;OJ] for rate rule [CO-CsH_O;O2b]
Euclidian distance = 1.0
Multiplied by reaction path degeneracy 2.0
family: R_Addition_MultipleBond
Ea raised from 119.3 to 123.5 kJ/mol to match endothermicity of reaction."""),
)

reaction(
    label = 'reaction7',
    reactants = ['[CH2]C[CH]C([O])O[O](8519)'],
    products = ['[CH2][CH]CC([O])O[O](8173)'],
    transitionState = 'TS7',
    kinetics = Arrhenius(A=(1.91383e+06,'s^-1'), n=2.02407, Ea=(137.004,'kJ/mol'), T0=(1,'K'), comment="""Estimated using an average for rate rule [R2H_S;C_rad_out_H/NonDeC;XH_out]
Euclidian distance = 0
Multiplied by reaction path degeneracy 2.0
family: intra_H_migration"""),
)

reaction(
    label = 'reaction8',
    reactants = ['[CH2][CH]CC([O])O[O](8173)'],
    products = ['[CH2][CH]C[C](O)O[O](8520)'],
    transitionState = 'TS8',
    kinetics = Arrhenius(A=(1.43381e+07,'s^-1'), n=1.70481, Ea=(158.036,'kJ/mol'), T0=(1,'K'), comment="""Estimated using template [R2H_S;Y_rad_out;Cs_H_out_NonDe] for rate rule [R2H_S;O_rad_out;Cs_H_out_NDMustO]
Euclidian distance = 1.41421356237
family: intra_H_migration"""),
)

reaction(
    label = 'reaction9',
    reactants = ['[CH2][CH]C[C]([O])OO(1366)'],
    products = ['[CH2][CH]CC([O])O[O](8173)'],
    transitionState = 'TS9',
    kinetics = Arrhenius(A=(40813.3,'s^-1'), n=2.17068, Ea=(135.172,'kJ/mol'), T0=(1,'K'), comment="""Estimated using average of templates [R3H_SS;Y_rad_out;O_H_out] + [R3H_SS_O;Y_rad_out;XH_out] for rate rule [R3H_SS_O;Y_rad_out;O_H_out]
Euclidian distance = 1.0
family: intra_H_migration"""),
)

reaction(
    label = 'reaction10',
    reactants = ['[CH2]CC[C]([O])O[O](1364)'],
    products = ['[CH2][CH]CC([O])O[O](8173)'],
    transitionState = 'TS10',
    kinetics = Arrhenius(A=(571202,'s^-1'), n=2.01653, Ea=(116.136,'kJ/mol'), T0=(1,'K'), comment="""Estimated using an average for rate rule [R3H_SS_Cs;Y_rad_out;XH_out]
Euclidian distance = 0
Multiplied by reaction path degeneracy 2.0
family: intra_H_migration"""),
)

reaction(
    label = 'reaction11',
    reactants = ['[CH2][CH]CC([O])O[O](8173)'],
    products = ['[CH2][CH][CH]C(O)O[O](8521)'],
    transitionState = 'TS11',
    kinetics = Arrhenius(A=(223829,'s^-1'), n=2.27675, Ea=(75.2806,'kJ/mol'), T0=(1,'K'), comment="""Estimated using an average for rate rule [R3H_SS_Cs;O_rad_out;XH_out]
Euclidian distance = 0
Multiplied by reaction path degeneracy 2.0
family: intra_H_migration"""),
)

reaction(
    label = 'reaction12',
    reactants = ['[CH2][CH]CC([O])O[O](8173)'],
    products = ['C[CH][CH]C([O])O[O](8279)'],
    transitionState = 'TS12',
    kinetics = Arrhenius(A=(2.50921e+06,'s^-1'), n=1.8375, Ea=(144.331,'kJ/mol'), T0=(1,'K'), comment="""Estimated using template [R3Hall;C_rad_out_2H;Cs_H_out_H/NonDeC] for rate rule [R3HJ;C_rad_out_2H;Cs_H_out_H/NonDeC]
Euclidian distance = 1.0
Multiplied by reaction path degeneracy 2.0
family: intra_H_migration"""),
)

reaction(
    label = 'reaction13',
    reactants = ['C[CH]C[C]([O])O[O](8280)'],
    products = ['[CH2][CH]CC([O])O[O](8173)'],
    transitionState = 'TS13',
    kinetics = Arrhenius(A=(190215,'s^-1'), n=1.95446, Ea=(128.468,'kJ/mol'), T0=(1,'K'), comment="""Estimated using template [R4Hall;Y_rad_out;Cs_H_out_2H] for rate rule [R4HJ_2;Y_rad_out;Cs_H_out_2H]
Euclidian distance = 1.0
Multiplied by reaction path degeneracy 3.0
family: intra_H_migration"""),
)

reaction(
    label = 'reaction14',
    reactants = ['[CH2][CH][CH]C([O])OO(8522)'],
    products = ['[CH2][CH]CC([O])O[O](8173)'],
    transitionState = 'TS14',
    kinetics = Arrhenius(A=(0.00504539,'s^-1'), n=3.83, Ea=(83.8892,'kJ/mol'), T0=(1,'K'), comment="""Estimated using an average for rate rule [R4H_SSS;Y_rad_out;O_H_out]
Euclidian distance = 0
family: intra_H_migration"""),
)

reaction(
    label = 'reaction15',
    reactants = ['O2(2)', '[CH2][CH]C[CH][O](749)'],
    products = ['[CH2][CH]CC([O])O[O](8173)'],
    transitionState = 'TS15',
    kinetics = Arrhenius(A=(2.18266e+06,'m^3/(mol*s)'), n=0.193158, Ea=(0,'kJ/mol'), T0=(1,'K'), comment="""Estimated using an average for rate rule [Y_rad;O2_birad]
Euclidian distance = 0
Multiplied by reaction path degeneracy 2.0
family: R_Recombination
Ea raised from -25.1 to 0 kJ/mol."""),
)

reaction(
    label = 'reaction16',
    reactants = ['[O][CH]O[O](8201)', '[CH2][CH][CH2](6136)'],
    products = ['[CH2][CH]CC([O])O[O](8173)'],
    transitionState = 'TS16',
    kinetics = Arrhenius(A=(3.9578e+07,'m^3/(mol*s)'), n=-0.126319, Ea=(0,'kJ/mol'), T0=(1,'K'), comment="""Estimated using an average for rate rule [Y_rad;Y_rad]
Euclidian distance = 0
Multiplied by reaction path degeneracy 2.0
family: R_Recombination
Ea raised from -15.6 to -15.6 kJ/mol.
Ea raised from -15.6 to 0 kJ/mol."""),
)

reaction(
    label = 'reaction17',
    reactants = ['H(8)', '[CH2][CH][CH]C([O])O[O](8523)'],
    products = ['[CH2][CH]CC([O])O[O](8173)'],
    transitionState = 'TS17',
    kinetics = Arrhenius(A=(4.34078e+06,'m^3/(mol*s)'), n=0.278577, Ea=(0,'kJ/mol'), T0=(1,'K'), comment="""Estimated using an average for rate rule [Y_rad;H_rad]
Euclidian distance = 0
family: R_Recombination
Ea raised from -1.4 to 0 kJ/mol."""),
)

reaction(
    label = 'reaction18',
    reactants = ['H(8)', '[CH2][CH]C[C]([O])O[O](8524)'],
    products = ['[CH2][CH]CC([O])O[O](8173)'],
    transitionState = 'TS18',
    kinetics = Arrhenius(A=(4.34078e+06,'m^3/(mol*s)'), n=0.278577, Ea=(0,'kJ/mol'), T0=(1,'K'), comment="""Estimated using an average for rate rule [Y_rad;H_rad]
Euclidian distance = 0
family: R_Recombination
Ea raised from -1.4 to 0 kJ/mol."""),
)

reaction(
    label = 'reaction19',
    reactants = ['[CH2][CH]CC([O])O[O](8173)'],
    products = ['C=C[CH]C(O)O[O](8525)'],
    transitionState = 'TS19',
    kinetics = Arrhenius(A=(1.4874e+09,'s^-1'), n=1.045, Ea=(63.4002,'kJ/mol'), T0=(1,'K'), comment="""Estimated using an average for rate rule [R3radExo;Y_rad;XH_Rrad]
Euclidian distance = 0
Multiplied by reaction path degeneracy 2.0
family: Intra_Disproportionation"""),
)

reaction(
    label = 'reaction20',
    reactants = ['[CH2][CH]CC([O])O[O](8173)'],
    products = ['C=CC[C]([O])OO(1285)'],
    transitionState = 'TS20',
    kinetics = Arrhenius(A=(7.437e+08,'s^-1'), n=1.045, Ea=(63.4002,'kJ/mol'), T0=(1,'K'), comment="""Estimated using an average for rate rule [R3radExo;Y_rad;XH_Rrad]
Euclidian distance = 0
family: Intra_Disproportionation"""),
)

reaction(
    label = 'reaction21',
    reactants = ['[CH2][CH]CC([O])O[O](8173)'],
    products = ['[CH2]CCC(=O)O[O](8526)'],
    transitionState = 'TS21',
    kinetics = Arrhenius(A=(7.437e+08,'s^-1'), n=1.045, Ea=(63.4002,'kJ/mol'), T0=(1,'K'), comment="""Estimated using an average for rate rule [R3radExo;Y_rad;XH_Rrad]
Euclidian distance = 0
family: Intra_Disproportionation"""),
)

reaction(
    label = 'reaction22',
    reactants = ['[CH2][CH]CC([O])O[O](8173)'],
    products = ['C=C[CH]C([O])OO(8527)'],
    transitionState = 'TS22',
    kinetics = Arrhenius(A=(1.02844e+09,'s^-1'), n=0.311, Ea=(24.9733,'kJ/mol'), T0=(1,'K'), comment="""Estimated using template [R4;Y_rad;XH_Rrad] for rate rule [R4radExo;Y_rad;XH_Rrad]
Euclidian distance = 1.0
Multiplied by reaction path degeneracy 2.0
family: Intra_Disproportionation"""),
)

reaction(
    label = 'reaction23',
    reactants = ['[CH2][CH]CC([O])O[O](8173)'],
    products = ['C=CC[CH][O](692)', 'O2(S)(5486)'],
    transitionState = 'TS23',
    kinetics = Arrhenius(A=(5e+12,'s^-1'), n=0, Ea=(0,'kJ/mol'), T0=(1,'K'), Tmin=(300,'K'), Tmax=(1500,'K'), comment="""Exact match found for rate rule [RJJ]
Euclidian distance = 0
family: 1,4_Linear_birad_scission"""),
)

reaction(
    label = 'reaction24',
    reactants = ['[CH2][CH]CC([O])O[O](8173)'],
    products = ['O(T)(63)', '[CH2][CH]CC1OO1(8528)'],
    transitionState = 'TS24',
    kinetics = Arrhenius(A=(2.30161e+09,'s^-1'), n=1.08, Ea=(189.53,'kJ/mol'), T0=(1,'K'), comment="""Estimated using an average for rate rule [R2OO_S;Y_rad_intra;OOJ]
Euclidian distance = 0
family: Cyclic_Ether_Formation"""),
)

reaction(
    label = 'reaction25',
    reactants = ['[CH2][CH]CC([O])O[O](8173)'],
    products = ['O(T)(63)', '[CH2]C1CC([O])O1(712)'],
    transitionState = 'TS25',
    kinetics = Arrhenius(A=(3.11355e+11,'s^-1'), n=0, Ea=(82.8981,'kJ/mol'), T0=(1,'K'), comment="""Estimated using template [R3OO_SS;Y_rad_intra;OO] for rate rule [R3OO_SS;Y_rad_intra;OOJ]
Euclidian distance = 1.0
family: Cyclic_Ether_Formation"""),
)

reaction(
    label = 'reaction26',
    reactants = ['[CH2]C([CH2])C([O])O[O](8175)'],
    products = ['[CH2][CH]CC([O])O[O](8173)'],
    transitionState = 'TS26',
    kinetics = Arrhenius(A=(2.66e+08,'s^-1'), n=1.36, Ea=(157.318,'kJ/mol'), T0=(1,'K'), comment="""Estimated using an average for rate rule [cCsCJ;CsJ-HH;C]
Euclidian distance = 0
Multiplied by reaction path degeneracy 2.0
family: 1,2_shiftC"""),
)

reaction(
    label = 'reaction27',
    reactants = ['[CH2][CH]CC([O])O[O](8173)'],
    products = ['[O]O(16)', '[CH2][CH]C=C[O](5307)'],
    transitionState = 'TS27',
    kinetics = Arrhenius(A=(1.91635e+11,'s^-1'), n=0.573333, Ea=(139.369,'kJ/mol'), T0=(1,'K'), comment="""Estimated using an average for rate rule [R2OO]
Euclidian distance = 0
Multiplied by reaction path degeneracy 2.0
family: HO2_Elimination_from_PeroxyRadical"""),
)

reaction(
    label = 'reaction28',
    reactants = ['[CH2][CH]CC([O])O[O](8173)'],
    products = ['[CH2][CH]CC1OOO1(8529)'],
    transitionState = 'TS28',
    kinetics = Arrhenius(A=(1.62e+12,'s^-1'), n=-0.305, Ea=(8.28432,'kJ/mol'), T0=(1,'K'), comment="""Estimated using template [R4_SSS;Y_rad_out;Ypri_rad_out] for rate rule [R4_SSS;O_rad;Opri_rad]
Euclidian distance = 1.41421356237
family: Birad_recombination"""),
)

reaction(
    label = 'reaction29',
    reactants = ['[CH2][CH]CC([O])O[O](8173)'],
    products = ['[CH2]C1CC(O[O])O1(8177)'],
    transitionState = 'TS29',
    kinetics = Arrhenius(A=(1.62e+12,'s^-1'), n=-0.305, Ea=(8.28432,'kJ/mol'), T0=(1,'K'), comment="""Estimated using template [R4_SSS;Y_rad_out;Ypri_rad_out] for rate rule [R4_SSS;Y_rad_out;Opri_rad]
Euclidian distance = 1.0
family: Birad_recombination"""),
)

reaction(
    label = 'reaction30',
    reactants = ['[CH2][CH]CC([O])O[O](8173)'],
    products = ['[CH2]C1CC([O])OO1(8530)'],
    transitionState = 'TS30',
    kinetics = Arrhenius(A=(7.76e+09,'s^-1'), n=0.311, Ea=(7.1128,'kJ/mol'), T0=(1,'K'), comment="""Estimated using template [R5_SSSS;Y_rad_out;Ypri_rad_out] for rate rule [R5_SSSS;O_rad;Ypri_rad_out]
Euclidian distance = 1.0
family: Birad_recombination"""),
)

reaction(
    label = 'reaction31',
    reactants = ['[CH2][CH]CC([O])O[O](8173)'],
    products = ['C=CCC([O])O[O](8172)'],
    transitionState = 'TS31',
    kinetics = Arrhenius(A=(6.31e+07,'s^-1'), n=0, Ea=(0,'kJ/mol'), T0=(1,'K'), Tmin=(300,'K'), Tmax=(1500,'K'), comment="""From training reaction 2 used for Y_12_10
Exact match found for rate rule [Y_12_10]
Euclidian distance = 0
family: 1,2-Birad_to_alkene"""),
)

reaction(
    label = 'reaction32',
    reactants = ['O(T)(63)', '[CH2][CH]CC([O])[O](2311)'],
    products = ['[CH2][CH]CC([O])O[O](8173)'],
    transitionState = 'TS32',
    kinetics = Arrhenius(A=(87.1677,'m^3/(mol*s)'), n=1.88017, Ea=(5.1666,'kJ/mol'), T0=(1,'K'), Tmin=(303.03,'K'), Tmax=(2000,'K'), comment="""From training reaction 3 used for O_rad/NonDe;O_birad
Exact match found for rate rule [O_rad/NonDe;O_birad]
Euclidian distance = 0
Multiplied by reaction path degeneracy 2.0
family: Birad_R_Recombination"""),
)

reaction(
    label = 'reaction33',
    reactants = ['[CH2]C([O])O[O](607)', '[CH][CH2](721)'],
    products = ['[CH2][CH]CC([O])O[O](8173)'],
    transitionState = 'TS33',
    kinetics = Arrhenius(A=(1.14854e+06,'m^3/(mol*s)'), n=0.575199, Ea=(34.3157,'kJ/mol'), T0=(1,'K'), comment="""Estimated using template [Y_rad;Birad] for rate rule [C_rad/H2/Cs;Birad]
Euclidian distance = 3.0
family: Birad_R_Recombination"""),
)

reaction(
    label = 'reaction34',
    reactants = ['O(T)(63)', '[CH2][CH]C[CH]O[O](8531)'],
    products = ['[CH2][CH]CC([O])O[O](8173)'],
    transitionState = 'TS34',
    kinetics = Arrhenius(A=(93609.6,'m^3/(mol*s)'), n=1.13083, Ea=(163.847,'kJ/mol'), T0=(1,'K'), Tmin=(303.03,'K'), Tmax=(2000,'K'), comment="""Estimated using template [Y_rad;O_birad] for rate rule [C_rad/H/CsO;O_birad]
Euclidian distance = 4.0
family: Birad_R_Recombination"""),
)

reaction(
    label = 'reaction35',
    reactants = ['H(8)', '[CH2][C]CC([O])O[O](8532)'],
    products = ['[CH2][CH]CC([O])O[O](8173)'],
    transitionState = 'TS35',
    kinetics = Arrhenius(A=(1e+07,'m^3/(mol*s)'), n=0, Ea=(0,'kJ/mol'), T0=(1,'K'), comment="""Estimated using an average for rate rule [H_rad;Birad]
Euclidian distance = 0
family: Birad_R_Recombination"""),
)

reaction(
    label = 'reaction36',
    reactants = ['H(8)', '[CH][CH]CC([O])O[O](8533)'],
    products = ['[CH2][CH]CC([O])O[O](8173)'],
    transitionState = 'TS36',
    kinetics = Arrhenius(A=(1e+07,'m^3/(mol*s)'), n=0, Ea=(0,'kJ/mol'), T0=(1,'K'), comment="""Estimated using an average for rate rule [H_rad;Birad]
Euclidian distance = 0
family: Birad_R_Recombination"""),
)

network(
    label = '2287',
    isomers = [
        '[CH2][CH]CC([O])O[O](8173)',
    ],
    reactants = [
        ('[CH2]C=C(87)', '[O]OC=O(5472)'),
    ],
    bathGas = {
        'N2': 0.25,
        'Ne': 0.25,
        'He': 0.25,
        'Ar': 0.25,
    },
)

pressureDependence(
    label = '2287',
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

