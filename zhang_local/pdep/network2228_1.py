species(
    label = '[O]OC([O])C[C]([O])[O](8116)',
    structure = SMILES('[O]OC([O])C[C]([O])[O]'),
    E0 = (265.231,'kJ/mol'),
    modes = [
        HarmonicOscillator(frequencies=([492.5,1135,1000,360,370,350,1380,1390,370,380,2900,435,2750,2850,1437.5,1250,1305,750,350,180,1791.4,1791.5,1791.58,1791.72],'cm^-1')),
        HinderedRotor(inertia=(0.258859,'amu*angstrom^2'), symmetry=1, barrier=(5.95168,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.258886,'amu*angstrom^2'), symmetry=1, barrier=(5.9523,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.258962,'amu*angstrom^2'), symmetry=1, barrier=(5.95404,'kJ/mol'), semiclassical=False),
    ],
    spinMultiplicity = 6,
    opticalIsomers = 1,
    molecularWeight = (119.053,'amu'),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[1.03242,0.0869866,-0.000188439,2.03915e-07,-7.94202e-11,31985.8,33.8285], Tmin=(100,'K'), Tmax=(858.853,'K')), NASAPolynomial(coeffs=[-2.09763,0.0478302,-2.62041e-05,5.13621e-09,-3.53859e-13,34505.3,59.9904], Tmin=(858.853,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(265.231,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(245.277,'J/(mol*K)'), comment="""Thermo group additivity estimation: group(O2s-OsCs) + group(O2s-CsH) + group(O2s-CsH) + group(O2s-CsH) + group(O2s-OsH) + group(Cs-CsCsHH) + group(Cs-CsOsOsH) + group(Cs-CsOsOsH) + radical(Cs_P) + radical(CCOJ) + radical(CCOJ) + radical(ROOJ) + radical(CCOJ)"""),
)

species(
    label = 'C=C([O])[O](1172)',
    structure = SMILES('C=C([O])[O]'),
    E0 = (-40.8548,'kJ/mol'),
    modes = [
        HarmonicOscillator(frequencies=([2950,3100,1380,975,1025,1650,350,440,435,1725,219.703,219.72],'cm^-1')),
    ],
    spinMultiplicity = 3,
    opticalIsomers = 1,
    molecularWeight = (58.0361,'amu'),
    collisionModel = TransportData(shapeIndex=2, epsilon=(3790.78,'J/mol'), sigma=(6.02099,'angstroms'), dipoleMoment=(0,'C*m'), polarizability=(0,'angstroms^3'), rotrelaxcollnum=0, comment="""Epsilon & sigma estimated with Tc=592.11 K, Pc=39.41 bar (from Joback method)"""),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[2.85742,0.0201906,-6.55878e-06,-8.62627e-09,5.39429e-12,-4868.13,11.9276], Tmin=(100,'K'), Tmax=(974.512,'K')), NASAPolynomial(coeffs=[9.20068,0.00568109,-1.96826e-06,3.71357e-10,-2.78204e-14,-6651.8,-21.3196], Tmin=(974.512,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(-40.8548,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(133.032,'J/(mol*K)'), comment="""Thermo group additivity estimation: group(O2s-(Cds-Cd)H) + group(O2s-(Cds-Cd)H) + group(Cds-CdsCsCs) + group(Cds-CdsHH) + radical(C=COJ) + radical(C=COJ)"""),
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
    label = '[O]OC([O])C=C([O])[O](9124)',
    structure = SMILES('[O]OC([O])C=C([O])[O]'),
    E0 = (12.2759,'kJ/mol'),
    modes = [
        HarmonicOscillator(frequencies=([492.5,1135,1000,1380,1390,370,380,2900,435,350,440,435,1725,3010,987.5,1337.5,450,1655,180,180,180,1302.28],'cm^-1')),
        HinderedRotor(inertia=(0.246978,'amu*angstrom^2'), symmetry=1, barrier=(5.6785,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.243507,'amu*angstrom^2'), symmetry=1, barrier=(5.5987,'kJ/mol'), semiclassical=False),
    ],
    spinMultiplicity = 5,
    opticalIsomers = 1,
    molecularWeight = (118.045,'amu'),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[0.816558,0.0804344,-0.0001461,1.36795e-07,-4.88854e-11,1581.04,29.3004], Tmin=(100,'K'), Tmax=(841.843,'K')), NASAPolynomial(coeffs=[7.32009,0.0282143,-1.50678e-05,2.94749e-09,-2.03911e-13,1241.48,3.53103], Tmin=(841.843,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(12.2759,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(224.491,'J/(mol*K)'), comment="""Thermo group additivity estimation: group(O2s-OsCs) + group(O2s-CsH) + group(O2s-(Cds-Cd)H) + group(O2s-(Cds-Cd)H) + group(O2s-OsH) + group(Cs-CsOsOsH) + group(Cds-CdsCsH) + group(Cds-CdsCsCs) + radical(ROOJ) + radical(CCOJ) + radical(C=COJ) + radical(C=COJ)"""),
)

species(
    label = '[O]O[C]([O])CC([O])=O(3308)',
    structure = SMILES('[O]O[C]([O])CC([O])=O'),
    E0 = (51.0983,'kJ/mol'),
    modes = [
        HarmonicOscillator(frequencies=([492.5,1135,1000,360,370,350,2750,2850,1437.5,1250,1305,750,350,200,800,933.333,1066.67,1200,1333.33,1466.67,1600],'cm^-1')),
        HinderedRotor(inertia=(0.156089,'amu*angstrom^2'), symmetry=1, barrier=(3.5888,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.156089,'amu*angstrom^2'), symmetry=1, barrier=(3.5888,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.156089,'amu*angstrom^2'), symmetry=1, barrier=(3.5888,'kJ/mol'), semiclassical=False),
    ],
    spinMultiplicity = 5,
    opticalIsomers = 1,
    molecularWeight = (118.045,'amu'),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[1.3351,0.0766949,-0.000160092,1.72771e-07,-6.77689e-11,6224.19,29.3586], Tmin=(100,'K'), Tmax=(850.702,'K')), NASAPolynomial(coeffs=[-1.18136,0.0440495,-2.41046e-05,4.74268e-09,-3.28427e-13,8261.76,50.5515], Tmin=(850.702,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(51.0983,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(220.334,'J/(mol*K)'), comment="""Thermo group additivity estimation: group(O2s-OsCs) + group(O2s-CsH) + group(O2s-(Cds-O2d)H) + group(O2s-OsH) + group(Cs-(Cds-O2d)CsHH) + group(Cs-CsOsOsH) + group(Cds-OdCsOs) + radical(CCOJ) + radical(Cs_P) + radical(CCOJ) + radical(ROOJ)"""),
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
    label = '[CH2][C]([O])[O](2304)',
    structure = SMILES('[CH2][C]([O])[O]'),
    E0 = (411.444,'kJ/mol'),
    modes = [
        HarmonicOscillator(frequencies=([360,370,350,3000,3100,440,815,1455,1000,1740.05,1740.26],'cm^-1')),
        HinderedRotor(inertia=(0.00349141,'amu*angstrom^2'), symmetry=1, barrier=(7.50228,'kJ/mol'), semiclassical=False),
    ],
    spinMultiplicity = 5,
    opticalIsomers = 1,
    molecularWeight = (58.0361,'amu'),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[3.01691,0.0315826,-7.13533e-05,8.28598e-08,-3.37644e-11,49511,16.8536], Tmin=(100,'K'), Tmax=(857.006,'K')), NASAPolynomial(coeffs=[-0.765291,0.0236786,-1.27872e-05,2.50409e-09,-1.7282e-13,51097.8,39.9925], Tmin=(857.006,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(411.444,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(128.874,'J/(mol*K)'), comment="""Thermo group additivity estimation: group(O2s-CsH) + group(O2s-CsH) + group(Cs-CsOsOsH) + group(Cs-CsHHH) + radical(CJCO) + radical(CCOJ) + radical(CCOJ) + radical(Cs_P)"""),
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
    label = '[O][CH]CC([O])=O(2298)',
    structure = SMILES('[O][CH]CC([O])=O'),
    E0 = (2.87112,'kJ/mol'),
    modes = [
        HarmonicOscillator(frequencies=([2750,2850,1437.5,1250,1305,750,350,3025,407.5,1350,352.5,180,180,1126.96,1126.98,1126.99,1127.06,1127.07,1127.17],'cm^-1')),
        HinderedRotor(inertia=(0.208927,'amu*angstrom^2'), symmetry=1, barrier=(4.80364,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.208846,'amu*angstrom^2'), symmetry=1, barrier=(4.80178,'kJ/mol'), semiclassical=False),
    ],
    spinMultiplicity = 4,
    opticalIsomers = 1,
    molecularWeight = (87.0541,'amu'),
    collisionModel = TransportData(shapeIndex=2, epsilon=(4425.74,'J/mol'), sigma=(6.62857,'angstroms'), dipoleMoment=(0,'C*m'), polarizability=(0,'angstroms^3'), rotrelaxcollnum=0, comment="""Epsilon & sigma estimated with Tc=691.29 K, Pc=34.48 bar (from Joback method)"""),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[2.13123,0.0562129,-0.000116907,1.31389e-07,-5.29461e-11,398.008,20.6838], Tmin=(100,'K'), Tmax=(853.341,'K')), NASAPolynomial(coeffs=[-2.56777,0.0400097,-2.12254e-05,4.13941e-09,-2.8566e-13,2591.89,50.7642], Tmin=(853.341,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(2.87112,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(199.547,'J/(mol*K)'), comment="""Thermo group additivity estimation: group(O2s-CsH) + group(O2s-(Cds-O2d)H) + group(Cs-(Cds-O2d)CsHH) + group(Cs-CsOsHH) + group(Cds-OdCsOs) + radical(CCOJ) + radical(CCsJOH) + radical(CCOJ)"""),
)

species(
    label = '[O]OC([O])[CH]C([O])[O](9125)',
    structure = SMILES('[O]OC([O])[CH]C([O])[O]'),
    E0 = (260.399,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
    molecularWeight = (119.053,'amu'),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[1.01827,0.0883042,-0.000193981,2.10938e-07,-8.22383e-11,31404.3,35.4963], Tmin=(100,'K'), Tmax=(860.654,'K')), NASAPolynomial(coeffs=[-2.62405,0.0488031,-2.67877e-05,5.24771e-09,-3.61151e-13,34121.2,64.6634], Tmin=(860.654,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(260.399,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(245.277,'J/(mol*K)'), comment="""Thermo group additivity estimation: group(O2s-OsCs) + group(O2s-CsH) + group(O2s-CsH) + group(O2s-CsH) + group(O2s-OsH) + group(Cs-CsCsHH) + group(Cs-CsOsOsH) + group(Cs-CsOsOsH) + radical(CCOJ) + radical(CCJCOOH) + radical(CCOJ) + radical(CCOJ) + radical(ROOJ)"""),
)

species(
    label = '[O]O[C](O)C[C]([O])[O](9126)',
    structure = SMILES('[O]O[C](O)C[C]([O])[O]'),
    E0 = (244.772,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
    molecularWeight = (119.053,'amu'),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[0.722936,0.0902608,-0.000185743,1.89727e-07,-7.10316e-11,29539.8,34.7681], Tmin=(100,'K'), Tmax=(862.069,'K')), NASAPolynomial(coeffs=[2.40133,0.0390949,-2.12367e-05,4.13882e-09,-2.83828e-13,30862.3,36.2687], Tmin=(862.069,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(244.772,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(241.12,'J/(mol*K)'), comment="""Thermo group additivity estimation: group(O2s-OsCs) + group(O2s-CsH) + group(O2s-CsH) + group(O2s-CsH) + group(O2s-OsH) + group(Cs-CsCsHH) + group(Cs-CsOsOsH) + group(Cs-CsOsOsH) + radical(Cs_P) + radical(ROOJ) + radical(Cs_P) + radical(CCOJ) + radical(CCOJ)"""),
)

species(
    label = '[O][C]([O])C[C]([O])OO(9127)',
    structure = SMILES('[O][C]([O])C[C]([O])OO'),
    E0 = (318.472,'kJ/mol'),
    modes = [
        HarmonicOscillator(frequencies=([3615,1310,387.5,850,1000,360,363.333,366.667,370,300,400,2750,2850,1437.5,1250,1305,750,350,215.676,810.486,1080.65,1385.59,1683.89],'cm^-1')),
        HinderedRotor(inertia=(0.130705,'amu*angstrom^2'), symmetry=1, barrier=(3.49471,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.130705,'amu*angstrom^2'), symmetry=1, barrier=(3.49471,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.130705,'amu*angstrom^2'), symmetry=1, barrier=(3.49471,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.130705,'amu*angstrom^2'), symmetry=1, barrier=(3.49471,'kJ/mol'), semiclassical=False),
    ],
    spinMultiplicity = 6,
    opticalIsomers = 1,
    molecularWeight = (119.053,'amu'),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[0.798967,0.0914074,-0.000194623,2.07654e-07,-8.05821e-11,38398.4,34.5946], Tmin=(100,'K'), Tmax=(850.893,'K')), NASAPolynomial(coeffs=[-0.641774,0.0470663,-2.63497e-05,5.21576e-09,-3.61814e-13,40493.9,52.1858], Tmin=(850.893,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(318.472,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(241.12,'J/(mol*K)'), comment="""Thermo group additivity estimation: group(O2s-OsCs) + group(O2s-CsH) + group(O2s-CsH) + group(O2s-CsH) + group(O2s-OsH) + group(Cs-CsCsHH) + group(Cs-CsOsOsH) + group(Cs-CsOsOsH) + radical(Cs_P) + radical(CCOJ) + radical(CCOJ) + radical(Cs_P) + radical(CCOJ)"""),
)

species(
    label = '[O]O[C]([O])CC([O])[O](9128)',
    structure = SMILES('[O]O[C]([O])CC([O])[O]'),
    E0 = (265.231,'kJ/mol'),
    modes = [
        HarmonicOscillator(frequencies=([492.5,1135,1000,360,370,350,1380,1390,370,380,2900,435,2750,2850,1437.5,1250,1305,750,350,180,1790.28,1791.57,1791.73,1792.76],'cm^-1')),
        HinderedRotor(inertia=(0.258432,'amu*angstrom^2'), symmetry=1, barrier=(5.94185,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.25887,'amu*angstrom^2'), symmetry=1, barrier=(5.95193,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.259385,'amu*angstrom^2'), symmetry=1, barrier=(5.96378,'kJ/mol'), semiclassical=False),
    ],
    spinMultiplicity = 6,
    opticalIsomers = 1,
    molecularWeight = (119.053,'amu'),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[1.03245,0.0869862,-0.000188437,2.03912e-07,-7.94192e-11,31985.8,33.8284], Tmin=(100,'K'), Tmax=(858.856,'K')), NASAPolynomial(coeffs=[-2.0977,0.0478304,-2.62042e-05,5.13622e-09,-3.5386e-13,34505.3,59.9907], Tmin=(858.856,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(265.231,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(245.277,'J/(mol*K)'), comment="""Thermo group additivity estimation: group(O2s-OsCs) + group(O2s-CsH) + group(O2s-CsH) + group(O2s-CsH) + group(O2s-OsH) + group(Cs-CsCsHH) + group(Cs-CsOsOsH) + group(Cs-CsOsOsH) + radical(CCOJ) + radical(ROOJ) + radical(CCOJ) + radical(Cs_P) + radical(CCOJ)"""),
)

species(
    label = '[O]OC(O)[CH][C]([O])[O](9129)',
    structure = SMILES('[O]OC(O)[CH][C]([O])[O]'),
    E0 = (239.94,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
    molecularWeight = (119.053,'amu'),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[0.708081,0.0915878,-0.000191324,1.9681e-07,-7.38804e-11,28958.3,36.4384], Tmin=(100,'K'), Tmax=(864.029,'K')), NASAPolynomial(coeffs=[1.87548,0.0400666,-2.18196e-05,4.25016e-09,-2.91106e-13,30478,40.9385], Tmin=(864.029,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(239.94,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(241.12,'J/(mol*K)'), comment="""Thermo group additivity estimation: group(O2s-OsCs) + group(O2s-CsH) + group(O2s-CsH) + group(O2s-CsH) + group(O2s-OsH) + group(Cs-CsCsHH) + group(Cs-CsOsOsH) + group(Cs-CsOsOsH) + radical(CCOJ) + radical(CCJCOOH) + radical(Cs_P) + radical(CCOJ) + radical(ROOJ)"""),
)

species(
    label = '[O]OC([O])[CH][C]([O])O(9130)',
    structure = SMILES('[O]OC([O])[CH][C]([O])O'),
    E0 = (239.94,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
    molecularWeight = (119.053,'amu'),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[0.708133,0.0915871,-0.000191322,1.96806e-07,-7.38785e-11,28958.3,37.1314], Tmin=(100,'K'), Tmax=(864.035,'K')), NASAPolynomial(coeffs=[1.87536,0.0400668,-2.18198e-05,4.25019e-09,-2.91109e-13,30478,41.6324], Tmin=(864.035,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(239.94,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(241.12,'J/(mol*K)'), comment="""Thermo group additivity estimation: group(O2s-OsCs) + group(O2s-CsH) + group(O2s-CsH) + group(O2s-CsH) + group(O2s-OsH) + group(Cs-CsCsHH) + group(Cs-CsOsOsH) + group(Cs-CsOsOsH) + radical(CCOJ) + radical(ROOJ) + radical(Cs_P) + radical(CCOJ) + radical(CCJCOOH)"""),
)

species(
    label = '[O]O[C]([O])C[C]([O])O(9131)',
    structure = SMILES('[O]O[C]([O])C[C]([O])O'),
    E0 = (244.772,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
    molecularWeight = (119.053,'amu'),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[0.723072,0.090259,-0.000185737,1.89717e-07,-7.10267e-11,29539.8,35.4608], Tmin=(100,'K'), Tmax=(862.085,'K')), NASAPolynomial(coeffs=[2.40103,0.0390954,-2.1237e-05,4.1389e-09,-2.83835e-13,30862.4,36.9635], Tmin=(862.085,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(244.772,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(241.12,'J/(mol*K)'), comment="""Thermo group additivity estimation: group(O2s-OsCs) + group(O2s-CsH) + group(O2s-CsH) + group(O2s-CsH) + group(O2s-OsH) + group(Cs-CsCsHH) + group(Cs-CsOsOsH) + group(Cs-CsOsOsH) + radical(ROOJ) + radical(CCOJ) + radical(Cs_P) + radical(Cs_P) + radical(CCOJ)"""),
)

species(
    label = '[O][C]([O])[CH]C([O])OO(9132)',
    structure = SMILES('[O][C]([O])[CH]C([O])OO'),
    E0 = (313.641,'kJ/mol'),
    modes = [
        HarmonicOscillator(frequencies=([3615,1310,387.5,850,1000,360,370,350,1380,1390,370,380,2900,435,3025,407.5,1350,352.5,244.611,855.935,1141.25,1482.04,1838.66],'cm^-1')),
        HinderedRotor(inertia=(0.0965648,'amu*angstrom^2'), symmetry=1, barrier=(3.32115,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.0965648,'amu*angstrom^2'), symmetry=1, barrier=(3.32115,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.0965648,'amu*angstrom^2'), symmetry=1, barrier=(3.32115,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.0965648,'amu*angstrom^2'), symmetry=1, barrier=(3.32115,'kJ/mol'), semiclassical=False),
    ],
    spinMultiplicity = 6,
    opticalIsomers = 1,
    molecularWeight = (119.053,'amu'),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[0.784379,0.0927299,-0.00020018,2.14689e-07,-8.34011e-11,37816.9,36.264], Tmin=(100,'K'), Tmax=(853.093,'K')), NASAPolynomial(coeffs=[-1.16456,0.0480328,-2.69296e-05,5.32638e-09,-3.69032e-13,40108.4,56.8385], Tmin=(853.093,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(313.641,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(241.12,'J/(mol*K)'), comment="""Thermo group additivity estimation: group(O2s-OsCs) + group(O2s-CsH) + group(O2s-CsH) + group(O2s-CsH) + group(O2s-OsH) + group(Cs-CsCsHH) + group(Cs-CsOsOsH) + group(Cs-CsOsOsH) + radical(CCOJ) + radical(CCOJ) + radical(Cs_P) + radical(CCOJ) + radical(CCJCOOH)"""),
)

species(
    label = '[O][CH]C[C]([O])[O](2619)',
    structure = SMILES('[O][CH]C[C]([O])[O]'),
    E0 = (422.25,'kJ/mol'),
    modes = [
        HarmonicOscillator(frequencies=([2750,2850,1437.5,1250,1305,750,350,360,370,350,3025,407.5,1350,352.5,180,1809.58,1809.7,1811.11,1811.71],'cm^-1')),
        HinderedRotor(inertia=(0.210269,'amu*angstrom^2'), symmetry=1, barrier=(4.8345,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.209899,'amu*angstrom^2'), symmetry=1, barrier=(4.82598,'kJ/mol'), semiclassical=False),
    ],
    spinMultiplicity = 6,
    opticalIsomers = 1,
    molecularWeight = (87.0541,'amu'),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[1.86238,0.0672735,-0.000155797,1.74833e-07,-6.90091e-11,50842.5,25.7204], Tmin=(100,'K'), Tmax=(867.888,'K')), NASAPolynomial(coeffs=[-3.65846,0.0414406,-2.25238e-05,4.3825e-09,-2.99825e-13,53732,62.6989], Tmin=(867.888,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(422.25,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(199.547,'J/(mol*K)'), comment="""Thermo group additivity estimation: group(O2s-CsH) + group(O2s-CsH) + group(O2s-CsH) + group(Cs-CsCsHH) + group(Cs-CsOsOsH) + group(Cs-CsOsHH) + radical(Cs_P) + radical(CCOJ) + radical(CCOJ) + radical(CCsJOH) + radical(CCOJ)"""),
)

species(
    label = '[O]OC([O])[CH][C]([O])[O](9133)',
    structure = SMILES('[O]OC([O])[CH][C]([O])[O]'),
    E0 = (465.646,'kJ/mol'),
    modes = [
        HarmonicOscillator(frequencies=([492.5,1135,1000,360,370,350,1380,1390,370,380,2900,435,3025,407.5,1350,352.5,180,1846.46,1846.56,1847.03,1847.07],'cm^-1')),
        HinderedRotor(inertia=(0.173396,'amu*angstrom^2'), symmetry=1, barrier=(3.98672,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.17346,'amu*angstrom^2'), symmetry=1, barrier=(3.9882,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.173185,'amu*angstrom^2'), symmetry=1, barrier=(3.98187,'kJ/mol'), semiclassical=False),
    ],
    spinMultiplicity = 7,
    opticalIsomers = 1,
    molecularWeight = (118.045,'amu'),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[1.05127,0.0890822,-0.00020455,2.23258e-07,-8.66493e-11,56087.2,36.066], Tmin=(100,'K'), Tmax=(865.543,'K')), NASAPolynomial(coeffs=[-2.79106,0.0464403,-2.5979e-05,5.09546e-09,-3.49733e-13,59014.7,67.1177], Tmin=(865.543,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(465.646,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(220.334,'J/(mol*K)'), comment="""Thermo group additivity estimation: group(O2s-OsCs) + group(O2s-CsH) + group(O2s-CsH) + group(O2s-CsH) + group(O2s-OsH) + group(Cs-CsCsHH) + group(Cs-CsOsOsH) + group(Cs-CsOsOsH) + radical(CCOJ) + radical(CCJCOOH) + radical(ROOJ) + radical(CCOJ) + radical(Cs_P) + radical(CCOJ)"""),
)

species(
    label = '[O]O[C]([O])C[C]([O])[O](9134)',
    structure = SMILES('[O]O[C]([O])C[C]([O])[O]'),
    E0 = (470.477,'kJ/mol'),
    modes = [
        HarmonicOscillator(frequencies=([492.5,1135,1000,360,363.333,366.667,370,300,400,2750,2850,1437.5,1250,1305,750,350,180,1865.91,1866.42,1866.77,1867.11],'cm^-1')),
        HinderedRotor(inertia=(0.19944,'amu*angstrom^2'), symmetry=1, barrier=(4.58552,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.199904,'amu*angstrom^2'), symmetry=1, barrier=(4.59618,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.199905,'amu*angstrom^2'), symmetry=1, barrier=(4.5962,'kJ/mol'), semiclassical=False),
    ],
    spinMultiplicity = 7,
    opticalIsomers = 1,
    molecularWeight = (118.045,'amu'),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[1.06539,0.0877646,-0.000199005,2.16226e-07,-8.38243e-11,56668.7,34.3984], Tmin=(100,'K'), Tmax=(864.067,'K')), NASAPolynomial(coeffs=[-2.26305,0.0454647,-2.53939e-05,4.98357e-09,-3.42409e-13,59398.2,62.4358], Tmin=(864.067,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(470.477,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(220.334,'J/(mol*K)'), comment="""Thermo group additivity estimation: group(O2s-OsCs) + group(O2s-CsH) + group(O2s-CsH) + group(O2s-CsH) + group(O2s-OsH) + group(Cs-CsCsHH) + group(Cs-CsOsOsH) + group(Cs-CsOsOsH) + radical(CCOJ) + radical(Cs_P) + radical(ROOJ) + radical(CCOJ) + radical(Cs_P) + radical(CCOJ)"""),
)

species(
    label = '[O]OC(O)C=C([O])[O](9135)',
    structure = SMILES('[O]OC(O)C=C([O])[O]'),
    E0 = (-213.429,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
    molecularWeight = (119.053,'amu'),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[0.734418,0.0792287,-0.000116143,8.11595e-08,-1.88913e-11,-25559,28.7722], Tmin=(100,'K'), Tmax=(633.137,'K')), NASAPolynomial(coeffs=[11.7478,0.0222876,-1.1184e-05,2.17052e-09,-1.51158e-13,-27206.9,-21.3277], Tmin=(633.137,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(-213.429,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(245.277,'J/(mol*K)'), comment="""Thermo group additivity estimation: group(O2s-OsCs) + group(O2s-CsH) + group(O2s-(Cds-Cd)H) + group(O2s-(Cds-Cd)H) + group(O2s-OsH) + group(Cs-CsOsOsH) + group(Cds-CdsCsH) + group(Cds-CdsCsCs) + radical(ROOJ) + radical(C=COJ) + radical(C=COJ)"""),
)

species(
    label = '[O][C](CC([O])=O)OO(3035)',
    structure = SMILES('[O][C](CC([O])=O)OO'),
    E0 = (-100.906,'kJ/mol'),
    modes = [
        HarmonicOscillator(frequencies=([3615,1310,387.5,850,1000,360,370,350,2750,2850,1437.5,1250,1305,750,350,200,800,933.333,1066.67,1200,1333.33,1466.67,1600],'cm^-1')),
        HinderedRotor(inertia=(0.156089,'amu*angstrom^2'), symmetry=1, barrier=(3.5888,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.156089,'amu*angstrom^2'), symmetry=1, barrier=(3.5888,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.156089,'amu*angstrom^2'), symmetry=1, barrier=(3.5888,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.156089,'amu*angstrom^2'), symmetry=1, barrier=(3.5888,'kJ/mol'), semiclassical=False),
    ],
    spinMultiplicity = 4,
    opticalIsomers = 1,
    molecularWeight = (119.053,'amu'),
    collisionModel = TransportData(shapeIndex=2, epsilon=(5118.62,'J/mol'), sigma=(7.52131,'angstroms'), dipoleMoment=(0,'C*m'), polarizability=(0,'angstroms^3'), rotrelaxcollnum=0, comment="""Epsilon & sigma estimated with Tc=799.52 K, Pc=27.3 bar (from Joback method)"""),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[1.07439,0.0802751,-0.000155529,1.6406e-07,-6.45341e-11,-12046.4,29.5341], Tmin=(100,'K'), Tmax=(832.276,'K')), NASAPolynomial(coeffs=[0.390336,0.0457378,-2.51115e-05,4.98711e-09,-3.48859e-13,-10622.5,40.5789], Tmin=(832.276,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(-100.906,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(241.12,'J/(mol*K)'), comment="""Thermo group additivity estimation: group(O2s-OsCs) + group(O2s-CsH) + group(O2s-(Cds-O2d)H) + group(O2s-OsH) + group(Cs-(Cds-O2d)CsHH) + group(Cs-CsOsOsH) + group(Cds-OdCsOs) + radical(CCOJ) + radical(Cs_P) + radical(CCOJ)"""),
)

species(
    label = '[O]OC(=O)CC([O])[O](9136)',
    structure = SMILES('[O]OC(=O)CC([O])[O]'),
    E0 = (-185.43,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
    molecularWeight = (119.053,'amu'),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[1.34354,0.0666283,-0.000104957,1.00139e-07,-3.81105e-11,-22214.4,30.1111], Tmin=(100,'K'), Tmax=(799.343,'K')), NASAPolynomial(coeffs=[4.24076,0.0364229,-1.87997e-05,3.69905e-09,-2.59708e-13,-22175.8,19.9215], Tmin=(799.343,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(-185.43,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(245.277,'J/(mol*K)'), comment="""Thermo group additivity estimation: group(O2s-CsH) + group(O2s-CsH) + group(O2s-O2s(Cds-O2d)) + group(O2s-OsH) + group(Cs-(Cds-O2d)CsHH) + group(Cs-CsOsOsH) + group(Cds-OdCsOs) + radical(CCOJ) + radical(C(=O)OOJ) + radical(CCOJ)"""),
)

species(
    label = '[O]C([O])=CC([O])OO(3302)',
    structure = SMILES('[O]C([O])=CC([O])OO'),
    E0 = (-139.729,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
    molecularWeight = (119.053,'amu'),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[0.579973,0.0837134,-0.000140416,1.26568e-07,-4.49987e-11,-16690.5,29.3906], Tmin=(100,'K'), Tmax=(797.949,'K')), NASAPolynomial(coeffs=[8.79761,0.0300705,-1.61748e-05,3.21616e-09,-2.26391e-13,-17605.7,-5.91651], Tmin=(797.949,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(-139.729,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(245.277,'J/(mol*K)'), comment="""Thermo group additivity estimation: group(O2s-OsCs) + group(O2s-CsH) + group(O2s-(Cds-Cd)H) + group(O2s-(Cds-Cd)H) + group(O2s-OsH) + group(Cs-CsOsOsH) + group(Cds-CdsCsH) + group(Cds-CdsCsCs) + radical(C=COJ) + radical(C=COJ) + radical(CCOJ)"""),
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
    label = '[O][C]([O])CC1OO1(9137)',
    structure = SMILES('[O][C]([O])CC1OO1'),
    E0 = (210.227,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
    molecularWeight = (103.054,'amu'),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[1.94845,0.0483926,-5.34835e-05,3.55092e-08,-1.00661e-11,25355.4,25.725], Tmin=(100,'K'), Tmax=(839.898,'K')), NASAPolynomial(coeffs=[6.87988,0.0249061,-1.15373e-05,2.21372e-09,-1.55322e-13,24527.1,2.79418], Tmin=(839.898,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(210.227,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(224.491,'J/(mol*K)'), comment="""Thermo group additivity estimation: group(O2s-OsCs) + group(O2s-OsCs) + group(O2s-CsH) + group(O2s-CsH) + group(Cs-CsCsHH) + group(Cs-CsOsOsH) + group(Cs-CsOsOsH) + ring(dioxirane) + radical(CCOJ) + radical(Cs_P) + radical(CCOJ)"""),
)

species(
    label = '[O]C1CC([O])([O])O1(9138)',
    structure = SMILES('[O]C1CC([O])([O])O1'),
    E0 = (-0.0440022,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
    molecularWeight = (103.054,'amu'),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[2.47801,0.0419322,-6.19732e-05,6.38131e-08,-2.43863e-11,41.5025,23.0007], Tmin=(100,'K'), Tmax=(905.247,'K')), NASAPolynomial(coeffs=[-0.912812,0.0352058,-1.4855e-05,2.62106e-09,-1.70795e-13,1544.92,43.935], Tmin=(905.247,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(-0.0440022,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(232.805,'J/(mol*K)'), comment="""Thermo group additivity estimation: group(O2s-CsCs) + group(O2s-CsH) + group(O2s-CsH) + group(O2s-CsH) + group(Cs-CsOsOsOs) + group(Cs-CsCsHH) + group(Cs-CsOsOsH) + ring(Oxetane) + radical(CCOJ) + radical(CCOJ) + radical(CCOJ)"""),
)

species(
    label = '[CH2]C([O])([O])C([O])O[O](8114)',
    structure = SMILES('[CH2]C([O])([O])C([O])O[O]'),
    E0 = (264.358,'kJ/mol'),
    modes = [
        HarmonicOscillator(frequencies=([492.5,1135,1000,1380,1390,370,380,2900,435,3000,3100,440,815,1455,1000,200,800,914.286,1028.57,1142.86,1257.14,1371.43,1485.71,1600],'cm^-1')),
        HinderedRotor(inertia=(0.156089,'amu*angstrom^2'), symmetry=1, barrier=(3.5888,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.156089,'amu*angstrom^2'), symmetry=1, barrier=(3.5888,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.156089,'amu*angstrom^2'), symmetry=1, barrier=(3.5888,'kJ/mol'), semiclassical=False),
    ],
    spinMultiplicity = 6,
    opticalIsomers = 1,
    molecularWeight = (119.053,'amu'),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[0.732074,0.0935495,-0.000199613,2.11311e-07,-8.10385e-11,31891.8,31.6285], Tmin=(100,'K'), Tmax=(860.686,'K')), NASAPolynomial(coeffs=[-0.575487,0.046658,-2.55775e-05,5.00393e-09,-3.44044e-13,34078.7,49.1378], Tmin=(860.686,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(264.358,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(245.277,'J/(mol*K)'), comment="""Thermo group additivity estimation: group(O2s-OsCs) + group(O2s-CsH) + group(O2s-CsH) + group(O2s-CsH) + group(O2s-OsH) + group(Cs-CsCsOsOs) + group(Cs-CsOsOsH) + group(Cs-CsHHH) + radical(CCOJ) + radical(ROOJ) + radical(CJC(O)2C) + radical(CC(C)(O)OJ) + radical(CC(C)(O)OJ)"""),
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
    label = '[O][CH]C=C([O])[O](3234)',
    structure = SMILES('[O][CH]C=C([O])[O]'),
    E0 = (112.984,'kJ/mol'),
    modes = [
        HarmonicOscillator(frequencies=([3025,407.5,1350,352.5,3010,987.5,1337.5,450,1655,350,440,435,1725,559.116,559.482,560.443,560.644],'cm^-1')),
        HinderedRotor(inertia=(0.000539674,'amu*angstrom^2'), symmetry=1, barrier=(0.119627,'kJ/mol'), semiclassical=False),
    ],
    spinMultiplicity = 5,
    opticalIsomers = 1,
    molecularWeight = (86.0462,'amu'),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[2.18591,0.0369791,-3.12753e-05,1.29059e-08,-2.12474e-12,13656.4,21.3605], Tmin=(100,'K'), Tmax=(1439.67,'K')), NASAPolynomial(coeffs=[10.9735,0.0125632,-5.83606e-06,1.12571e-09,-7.90829e-14,11126.2,-24.237], Tmin=(1439.67,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(112.984,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(178.761,'J/(mol*K)'), comment="""Thermo group additivity estimation: group(O2s-CsH) + group(O2s-(Cds-Cd)H) + group(O2s-(Cds-Cd)H) + group(Cs-(Cds-Cds)OsHH) + group(Cds-CdsCsH) + group(Cds-CdsCsCs) + radical(C=COJ) + radical(C=COJ) + radical(CCOJ) + radical(C=CCJO)"""),
)

species(
    label = '[O][C]([O])CC1OOO1(9139)',
    structure = SMILES('[O][C]([O])CC1OOO1'),
    E0 = (251.439,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
    molecularWeight = (119.053,'amu'),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[1.70284,0.0457634,-2.76147e-05,6.42277e-09,-4.14858e-13,30327.8,29.7315], Tmin=(100,'K'), Tmax=(1680.4,'K')), NASAPolynomial(coeffs=[16.6721,0.017659,-9.24736e-06,1.80188e-09,-1.24027e-13,24234.1,-53.4184], Tmin=(1680.4,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(251.439,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(249.434,'J/(mol*K)'), comment="""Thermo group additivity estimation: group(O2s-OsCs) + group(O2s-OsCs) + group(O2s-CsH) + group(O2s-CsH) + group(O2s-OsOs) + group(Cs-CsCsHH) + group(Cs-CsOsOsH) + group(Cs-CsOsOsH) + ring(Cyclobutane) + radical(Cs_P) + radical(CCOJ) + radical(CCOJ)"""),
)

species(
    label = '[O]OC1CC([O])([O])O1(8118)',
    structure = SMILES('[O]OC1CC([O])([O])O1'),
    E0 = (-2.23946,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
    molecularWeight = (119.053,'amu'),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[1.68711,0.0582253,-8.63335e-05,7.75605e-08,-2.65881e-11,-192.974,26.4991], Tmin=(100,'K'), Tmax=(923.088,'K')), NASAPolynomial(coeffs=[3.35763,0.0319929,-1.28422e-05,2.19352e-09,-1.39556e-13,307.833,22.9567], Tmin=(923.088,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(-2.23946,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(253.591,'J/(mol*K)'), comment="""Thermo group additivity estimation: group(O2s-CsCs) + group(O2s-OsCs) + group(O2s-CsH) + group(O2s-CsH) + group(O2s-OsH) + group(Cs-CsOsOsOs) + group(Cs-CsCsHH) + group(Cs-CsOsOsH) + ring(Oxetane) + radical(ROOJ) + radical(CCOJ) + radical(CCOJ)"""),
)

species(
    label = '[O]C1CC([O])([O])OO1(9140)',
    structure = SMILES('[O]C1CC([O])([O])OO1'),
    E0 = (-28.145,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
    molecularWeight = (119.053,'amu'),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[1.91071,0.0514405,-6.62202e-05,5.81223e-08,-2.04094e-11,-3315.07,25.8235], Tmin=(100,'K'), Tmax=(891.626,'K')), NASAPolynomial(coeffs=[2.97531,0.0335062,-1.39123e-05,2.46021e-09,-1.6172e-13,-2981.88,23.7426], Tmin=(891.626,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(-28.145,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(257.749,'J/(mol*K)'), comment="""Thermo group additivity estimation: group(O2s-OsCs) + group(O2s-OsCs) + group(O2s-CsH) + group(O2s-CsH) + group(O2s-CsH) + group(Cs-CsCsHH) + group(Cs-CsOsOsOs) + group(Cs-CsOsOsH) + ring(12dioxolane) + radical(CCOJ) + radical(CCOJ) + radical(CCOJ)"""),
)

species(
    label = '[O]OC([O])CC([O])=O(3303)',
    structure = SMILES('[O]OC([O])CC([O])=O'),
    E0 = (-154.148,'kJ/mol'),
    modes = [
        HarmonicOscillator(frequencies=([492.5,1135,1000,1380,1390,370,380,2900,435,2750,2850,1437.5,1250,1305,750,350,200,800,933.333,1066.67,1200,1333.33,1466.67,1600],'cm^-1')),
        HinderedRotor(inertia=(0.156089,'amu*angstrom^2'), symmetry=1, barrier=(3.5888,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.156089,'amu*angstrom^2'), symmetry=1, barrier=(3.5888,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.156089,'amu*angstrom^2'), symmetry=1, barrier=(3.5888,'kJ/mol'), semiclassical=False),
    ],
    spinMultiplicity = 4,
    opticalIsomers = 1,
    molecularWeight = (119.053,'amu'),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[1.30492,0.0758851,-0.000149425,1.60359e-07,-6.33436e-11,-18458.8,28.7787], Tmin=(100,'K'), Tmax=(843.155,'K')), NASAPolynomial(coeffs=[-1.03623,0.0464505,-2.49357e-05,4.90032e-09,-3.40296e-13,-16622.9,48.2197], Tmin=(843.155,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(-154.148,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(245.277,'J/(mol*K)'), comment="""Thermo group additivity estimation: group(O2s-OsCs) + group(O2s-CsH) + group(O2s-(Cds-O2d)H) + group(O2s-OsH) + group(Cs-(Cds-O2d)CsHH) + group(Cs-CsOsOsH) + group(Cds-OdCsOs) + radical(ROOJ) + radical(CCOJ) + radical(CCOJ)"""),
)

species(
    label = '[O][C]([O])CC([O])[O](9141)',
    structure = SMILES('[O][C]([O])CC([O])[O]'),
    E0 = (267.426,'kJ/mol'),
    modes = [
        HarmonicOscillator(frequencies=([2750,2850,1437.5,1250,1305,750,350,1380,1390,370,380,2900,435,360,370,350,180,1856.02,1856.36,1856.41,1856.47,1856.88],'cm^-1')),
        HinderedRotor(inertia=(0.15512,'amu*angstrom^2'), symmetry=1, barrier=(3.56651,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.155173,'amu*angstrom^2'), symmetry=1, barrier=(3.56774,'kJ/mol'), semiclassical=False),
    ],
    spinMultiplicity = 6,
    opticalIsomers = 1,
    molecularWeight = (103.054,'amu'),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[1.80698,0.0709141,-0.000165004,1.91642e-07,-7.79988e-11,32221,29.6939], Tmin=(100,'K'), Tmax=(855.358,'K')), NASAPolynomial(coeffs=[-6.3822,0.0510646,-2.82282e-05,5.56625e-09,-3.85293e-13,35749,80.3563], Tmin=(855.358,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(267.426,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(224.491,'J/(mol*K)'), comment="""Thermo group additivity estimation: group(O2s-CsH) + group(O2s-CsH) + group(O2s-CsH) + group(O2s-CsH) + group(Cs-CsCsHH) + group(Cs-CsOsOsH) + group(Cs-CsOsOsH) + radical(Cs_P) + radical(CCOJ) + radical(CCOJ) + radical(CCOJ) + radical(CCOJ)"""),
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
    label = '[O][C][O](2319)',
    structure = SMILES('[O][C][O]'),
    E0 = (504.804,'kJ/mol'),
    modes = [
        HarmonicOscillator(frequencies=([786.798,787.166,787.343],'cm^-1')),
    ],
    spinMultiplicity = 5,
    opticalIsomers = 1,
    molecularWeight = (44.0095,'amu'),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[4.36312,0.000897068,5.83704e-06,-5.0699e-09,1.05455e-12,60693,5.58445], Tmin=(100,'K'), Tmax=(1870.77,'K')), NASAPolynomial(coeffs=[6.74008,0.00434757,-3.77128e-06,7.92205e-10,-5.4643e-14,58310.5,-11.3625], Tmin=(1870.77,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(504.804,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(58.2013,'J/(mol*K)'), comment="""Thermo group additivity estimation: group(O2s-CsH) + group(O2s-CsH) + group(Cs-OsOsHH) + radical(CH2_triplet) + radical(OCOJ) + radical(OCOJ)"""),
)

species(
    label = '[O]O[CH]C[C]([O])[O](9142)',
    structure = SMILES('[O]O[CH]C[C]([O])[O]'),
    E0 = (428.339,'kJ/mol'),
    modes = [
        HarmonicOscillator(frequencies=([3025,407.5,1350,352.5,492.5,1135,1000,360,370,350,2750,2850,1437.5,1250,1305,750,350,180,1638.73,1639.15,1639.71],'cm^-1')),
        HinderedRotor(inertia=(0.204992,'amu*angstrom^2'), symmetry=1, barrier=(4.71318,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.204606,'amu*angstrom^2'), symmetry=1, barrier=(4.7043,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.204474,'amu*angstrom^2'), symmetry=1, barrier=(4.70127,'kJ/mol'), semiclassical=False),
    ],
    spinMultiplicity = 6,
    opticalIsomers = 1,
    molecularWeight = (103.054,'amu'),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[1.31679,0.0778839,-0.000166853,1.776e-07,-6.80635e-11,51595.8,28.9392], Tmin=(100,'K'), Tmax=(867.841,'K')), NASAPolynomial(coeffs=[-0.677308,0.040477,-2.16562e-05,4.19089e-09,-2.86029e-13,53696.7,48.3868], Tmin=(867.841,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(428.339,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(220.334,'J/(mol*K)'), comment="""Thermo group additivity estimation: group(O2s-OsCs) + group(O2s-CsH) + group(O2s-CsH) + group(O2s-OsH) + group(Cs-CsCsHH) + group(Cs-CsOsOsH) + group(Cs-CsOsHH) + radical(CCOJ) + radical(Cs_P) + radical(CCOJ) + radical(CCsJOOH) + radical(ROOJ)"""),
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
    E0 = (265.231,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS2',
    E0 = (265.231,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS3',
    E0 = (309.408,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS4',
    E0 = (265.231,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS5',
    E0 = (273.103,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS6',
    E0 = (265.231,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS7',
    E0 = (422.104,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS8',
    E0 = (423.267,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS9',
    E0 = (453.644,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS10',
    E0 = (381.367,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS11',
    E0 = (340.511,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS12',
    E0 = (420.7,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS13',
    E0 = (388.507,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS14',
    E0 = (397.53,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS15',
    E0 = (413.623,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS16',
    E0 = (640.332,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS17',
    E0 = (677.45,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS18',
    E0 = (682.282,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS19',
    E0 = (328.631,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS20',
    E0 = (328.631,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS21',
    E0 = (328.631,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS22',
    E0 = (290.204,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS23',
    E0 = (265.231,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS24',
    E0 = (454.761,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS25',
    E0 = (348.129,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS26',
    E0 = (421.676,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS27',
    E0 = (404.6,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS28',
    E0 = (273.515,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS29',
    E0 = (273.515,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS30',
    E0 = (272.343,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS31',
    E0 = (265.231,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS32',
    E0 = (515.627,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS33',
    E0 = (745.495,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS34',
    E0 = (835.22,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

reaction(
    label = 'reaction1',
    reactants = ['[O]OC([O])C[C]([O])[O](8116)'],
    products = ['C=C([O])[O](1172)', '[O]OC=O(5472)'],
    transitionState = 'TS1',
    kinetics = Arrhenius(A=(5e+12,'s^-1'), n=0, Ea=(0,'kJ/mol'), T0=(1,'K'), Tmin=(300,'K'), Tmax=(1500,'K'), comment="""Exact match found for rate rule [RJJ]
Euclidian distance = 0
family: 1,4_Linear_birad_scission"""),
)

reaction(
    label = 'reaction2',
    reactants = ['H(8)', '[O]OC([O])C=C([O])[O](9124)'],
    products = ['[O]OC([O])C[C]([O])[O](8116)'],
    transitionState = 'TS2',
    kinetics = Arrhenius(A=(8.22e+08,'cm^3/(mol*s)'), n=1.533, Ea=(41.1501,'kJ/mol'), T0=(1,'K'), comment="""From training reaction 192 used for Cd_R;HJ
Exact match found for rate rule [Cd_R;HJ]
Euclidian distance = 0
family: R_Addition_MultipleBond
Ea raised from 35.6 to 41.2 kJ/mol to match endothermicity of reaction."""),
)

reaction(
    label = 'reaction3',
    reactants = ['H(8)', '[O]O[C]([O])CC([O])=O(3308)'],
    products = ['[O]OC([O])C[C]([O])[O](8116)'],
    transitionState = 'TS3',
    kinetics = Arrhenius(A=(1.83701,'m^3/(mol*s)'), n=1.71338, Ea=(46.5052,'kJ/mol'), T0=(1,'K'), comment="""Estimated using average of templates [CO_O;HJ] + [CO-NdNd_O;YJ] for rate rule [CO-NdNd_O;HJ]
Euclidian distance = 1.0
family: R_Addition_MultipleBond"""),
)

reaction(
    label = 'reaction4',
    reactants = ['C=C([O])[O](1172)', '[O][CH]O[O](8201)'],
    products = ['[O]OC([O])C[C]([O])[O](8116)'],
    transitionState = 'TS4',
    kinetics = Arrhenius(A=(0.123469,'m^3/(mol*s)'), n=2.00579, Ea=(77.1976,'kJ/mol'), T0=(1,'K'), comment="""Estimated using an average for rate rule [Cd_R;CJ]
Euclidian distance = 0
family: R_Addition_MultipleBond
Ea raised from 72.0 to 77.2 kJ/mol to match endothermicity of reaction."""),
)

reaction(
    label = 'reaction5',
    reactants = ['[CH2][C]([O])[O](2304)', '[O]OC=O(5472)'],
    products = ['[O]OC([O])C[C]([O])[O](8116)'],
    transitionState = 'TS5',
    kinetics = Arrhenius(A=(6.87291,'m^3/(mol*s)'), n=1.39198, Ea=(57.1551,'kJ/mol'), T0=(1,'K'), comment="""Estimated using average of templates [CO_O;CJ] + [CO-NdH_O;YJ] for rate rule [CO-NdH_O;CJ]
Euclidian distance = 1.0
family: R_Addition_MultipleBond"""),
)

reaction(
    label = 'reaction6',
    reactants = ['O2(2)', '[O][CH]CC([O])=O(2298)'],
    products = ['[O]OC([O])C[C]([O])[O](8116)'],
    transitionState = 'TS6',
    kinetics = Arrhenius(A=(0.0849,'cm^3/(mol*s)'), n=3.486, Ea=(270.986,'kJ/mol'), T0=(1,'K'), Tmin=(300,'K'), Tmax=(1600,'K'), comment="""Estimated using template [CO-CsH_O;OJ] for rate rule [CO-CsH_O;O2b]
Euclidian distance = 1.0
Multiplied by reaction path degeneracy 2.0
family: R_Addition_MultipleBond
Ea raised from 267.1 to 271.0 kJ/mol to match endothermicity of reaction."""),
)

reaction(
    label = 'reaction7',
    reactants = ['[O]OC([O])C[C]([O])[O](8116)'],
    products = ['[O]OC([O])[CH]C([O])[O](9125)'],
    transitionState = 'TS7',
    kinetics = Arrhenius(A=(791180,'s^-1'), n=2.19286, Ea=(156.873,'kJ/mol'), T0=(1,'K'), comment="""Estimated using an average for rate rule [R2H_S;Y_rad_out;Cs_H_out_H/NonDeC]
Euclidian distance = 0
Multiplied by reaction path degeneracy 2.0
family: intra_H_migration"""),
)

reaction(
    label = 'reaction8',
    reactants = ['[O]OC([O])C[C]([O])[O](8116)'],
    products = ['[O]O[C](O)C[C]([O])[O](9126)'],
    transitionState = 'TS8',
    kinetics = Arrhenius(A=(1.43381e+07,'s^-1'), n=1.70481, Ea=(158.036,'kJ/mol'), T0=(1,'K'), comment="""Estimated using template [R2H_S;Y_rad_out;Cs_H_out_NonDe] for rate rule [R2H_S;O_rad_out;Cs_H_out_NDMustO]
Euclidian distance = 1.41421356237
family: intra_H_migration"""),
)

reaction(
    label = 'reaction9',
    reactants = ['[O][C]([O])C[C]([O])OO(9127)'],
    products = ['[O]OC([O])C[C]([O])[O](8116)'],
    transitionState = 'TS9',
    kinetics = Arrhenius(A=(40813.3,'s^-1'), n=2.17068, Ea=(135.172,'kJ/mol'), T0=(1,'K'), comment="""Estimated using average of templates [R3H_SS;Y_rad_out;O_H_out] + [R3H_SS_O;Y_rad_out;XH_out] for rate rule [R3H_SS_O;Y_rad_out;O_H_out]
Euclidian distance = 1.0
family: intra_H_migration"""),
)

reaction(
    label = 'reaction10',
    reactants = ['[O]O[C]([O])CC([O])[O](9128)'],
    products = ['[O]OC([O])C[C]([O])[O](8116)'],
    transitionState = 'TS10',
    kinetics = Arrhenius(A=(285601,'s^-1'), n=2.01653, Ea=(116.136,'kJ/mol'), T0=(1,'K'), comment="""Estimated using an average for rate rule [R3H_SS_Cs;Y_rad_out;XH_out]
Euclidian distance = 0
family: intra_H_migration"""),
)

reaction(
    label = 'reaction11',
    reactants = ['[O]OC([O])C[C]([O])[O](8116)'],
    products = ['[O]OC(O)[CH][C]([O])[O](9129)'],
    transitionState = 'TS11',
    kinetics = Arrhenius(A=(223829,'s^-1'), n=2.27675, Ea=(75.2806,'kJ/mol'), T0=(1,'K'), comment="""Estimated using an average for rate rule [R3H_SS_Cs;O_rad_out;XH_out]
Euclidian distance = 0
Multiplied by reaction path degeneracy 2.0
family: intra_H_migration"""),
)

reaction(
    label = 'reaction12',
    reactants = ['[O]OC([O])C[C]([O])[O](8116)'],
    products = ['[O]OC([O])[CH][C]([O])O(9130)'],
    transitionState = 'TS12',
    kinetics = Arrhenius(A=(2.61628e+09,'s^-1'), n=1.19923, Ea=(155.469,'kJ/mol'), T0=(1,'K'), comment="""Estimated using template [R3Hall;O_rad_out;Cs_H_out_H/NonDeC] for rate rule [R3HJ;O_rad_out;Cs_H_out_H/NonDeC]
Euclidian distance = 1.0
Multiplied by reaction path degeneracy 4.0
family: intra_H_migration"""),
)

reaction(
    label = 'reaction13',
    reactants = ['[O]OC([O])C[C]([O])[O](8116)'],
    products = ['[O]O[C]([O])C[C]([O])O(9131)'],
    transitionState = 'TS13',
    kinetics = Arrhenius(A=(1.04154e+06,'s^-1'), n=1.9431, Ea=(123.276,'kJ/mol'), T0=(1,'K'), comment="""Estimated using template [R4Hall;O_rad_out;XH_out] for rate rule [R4HJ_1;O_rad_out;XH_out]
Euclidian distance = 1.0
Multiplied by reaction path degeneracy 2.0
family: intra_H_migration"""),
)

reaction(
    label = 'reaction14',
    reactants = ['[O][C]([O])[CH]C([O])OO(9132)'],
    products = ['[O]OC([O])C[C]([O])[O](8116)'],
    transitionState = 'TS14',
    kinetics = Arrhenius(A=(0.00504539,'s^-1'), n=3.83, Ea=(83.8892,'kJ/mol'), T0=(1,'K'), comment="""Estimated using an average for rate rule [R4H_SSS;Y_rad_out;O_H_out]
Euclidian distance = 0
family: intra_H_migration"""),
)

reaction(
    label = 'reaction15',
    reactants = ['O2(2)', '[O][CH]C[C]([O])[O](2619)'],
    products = ['[O]OC([O])C[C]([O])[O](8116)'],
    transitionState = 'TS15',
    kinetics = Arrhenius(A=(2.18266e+06,'m^3/(mol*s)'), n=0.193158, Ea=(0,'kJ/mol'), T0=(1,'K'), comment="""Estimated using an average for rate rule [Y_rad;O2_birad]
Euclidian distance = 0
Multiplied by reaction path degeneracy 2.0
family: R_Recombination
Ea raised from -25.1 to 0 kJ/mol."""),
)

reaction(
    label = 'reaction16',
    reactants = ['[CH2][C]([O])[O](2304)', '[O][CH]O[O](8201)'],
    products = ['[O]OC([O])C[C]([O])[O](8116)'],
    transitionState = 'TS16',
    kinetics = Arrhenius(A=(1.9789e+07,'m^3/(mol*s)'), n=-0.126319, Ea=(0,'kJ/mol'), T0=(1,'K'), comment="""Estimated using an average for rate rule [Y_rad;Y_rad]
Euclidian distance = 0
family: R_Recombination
Ea raised from -15.6 to -15.6 kJ/mol.
Ea raised from -15.6 to 0 kJ/mol."""),
)

reaction(
    label = 'reaction17',
    reactants = ['H(8)', '[O]OC([O])[CH][C]([O])[O](9133)'],
    products = ['[O]OC([O])C[C]([O])[O](8116)'],
    transitionState = 'TS17',
    kinetics = Arrhenius(A=(4.34078e+06,'m^3/(mol*s)'), n=0.278577, Ea=(0,'kJ/mol'), T0=(1,'K'), comment="""Estimated using an average for rate rule [Y_rad;H_rad]
Euclidian distance = 0
family: R_Recombination
Ea raised from -1.4 to 0 kJ/mol."""),
)

reaction(
    label = 'reaction18',
    reactants = ['H(8)', '[O]O[C]([O])C[C]([O])[O](9134)'],
    products = ['[O]OC([O])C[C]([O])[O](8116)'],
    transitionState = 'TS18',
    kinetics = Arrhenius(A=(4.34078e+06,'m^3/(mol*s)'), n=0.278577, Ea=(0,'kJ/mol'), T0=(1,'K'), comment="""Estimated using an average for rate rule [H_rad;Y_rad]
Euclidian distance = 0
family: R_Recombination
Ea raised from -1.4 to 0 kJ/mol."""),
)

reaction(
    label = 'reaction19',
    reactants = ['[O]OC([O])C[C]([O])[O](8116)'],
    products = ['[O]OC(O)C=C([O])[O](9135)'],
    transitionState = 'TS19',
    kinetics = Arrhenius(A=(1.4874e+09,'s^-1'), n=1.045, Ea=(63.4002,'kJ/mol'), T0=(1,'K'), comment="""Estimated using an average for rate rule [R3radExo;Y_rad;XH_Rrad]
Euclidian distance = 0
Multiplied by reaction path degeneracy 2.0
family: Intra_Disproportionation"""),
)

reaction(
    label = 'reaction20',
    reactants = ['[O]OC([O])C[C]([O])[O](8116)'],
    products = ['[O][C](CC([O])=O)OO(3035)'],
    transitionState = 'TS20',
    kinetics = Arrhenius(A=(7.437e+08,'s^-1'), n=1.045, Ea=(63.4002,'kJ/mol'), T0=(1,'K'), comment="""Estimated using an average for rate rule [R3radExo;Y_rad;XH_Rrad]
Euclidian distance = 0
family: Intra_Disproportionation"""),
)

reaction(
    label = 'reaction21',
    reactants = ['[O]OC([O])C[C]([O])[O](8116)'],
    products = ['[O]OC(=O)CC([O])[O](9136)'],
    transitionState = 'TS21',
    kinetics = Arrhenius(A=(7.437e+08,'s^-1'), n=1.045, Ea=(63.4002,'kJ/mol'), T0=(1,'K'), comment="""Estimated using an average for rate rule [R3radExo;Y_rad;XH_Rrad]
Euclidian distance = 0
family: Intra_Disproportionation"""),
)

reaction(
    label = 'reaction22',
    reactants = ['[O]OC([O])C[C]([O])[O](8116)'],
    products = ['[O]C([O])=CC([O])OO(3302)'],
    transitionState = 'TS22',
    kinetics = Arrhenius(A=(1.02844e+09,'s^-1'), n=0.311, Ea=(24.9733,'kJ/mol'), T0=(1,'K'), comment="""Estimated using template [R4;Y_rad;XH_Rrad] for rate rule [R4radExo;Y_rad;XH_Rrad]
Euclidian distance = 1.0
Multiplied by reaction path degeneracy 2.0
family: Intra_Disproportionation"""),
)

reaction(
    label = 'reaction23',
    reactants = ['[O]OC([O])C[C]([O])[O](8116)'],
    products = ['[O][CH]CC([O])=O(2298)', 'O2(S)(5486)'],
    transitionState = 'TS23',
    kinetics = Arrhenius(A=(5e+12,'s^-1'), n=0, Ea=(0,'kJ/mol'), T0=(1,'K'), Tmin=(300,'K'), Tmax=(1500,'K'), comment="""Exact match found for rate rule [RJJ]
Euclidian distance = 0
family: 1,4_Linear_birad_scission"""),
)

reaction(
    label = 'reaction24',
    reactants = ['[O]OC([O])C[C]([O])[O](8116)'],
    products = ['O(T)(63)', '[O][C]([O])CC1OO1(9137)'],
    transitionState = 'TS24',
    kinetics = Arrhenius(A=(2.30161e+09,'s^-1'), n=1.08, Ea=(189.53,'kJ/mol'), T0=(1,'K'), comment="""Estimated using an average for rate rule [R2OO_S;Y_rad_intra;OOJ]
Euclidian distance = 0
family: Cyclic_Ether_Formation"""),
)

reaction(
    label = 'reaction25',
    reactants = ['[O]OC([O])C[C]([O])[O](8116)'],
    products = ['O(T)(63)', '[O]C1CC([O])([O])O1(9138)'],
    transitionState = 'TS25',
    kinetics = Arrhenius(A=(3.11355e+11,'s^-1'), n=0, Ea=(82.8981,'kJ/mol'), T0=(1,'K'), comment="""Estimated using template [R3OO_SS;Y_rad_intra;OO] for rate rule [R3OO_SS;Y_rad_intra;OOJ]
Euclidian distance = 1.0
family: Cyclic_Ether_Formation"""),
)

reaction(
    label = 'reaction26',
    reactants = ['[CH2]C([O])([O])C([O])O[O](8114)'],
    products = ['[O]OC([O])C[C]([O])[O](8116)'],
    transitionState = 'TS26',
    kinetics = Arrhenius(A=(1.33e+08,'s^-1'), n=1.36, Ea=(157.318,'kJ/mol'), T0=(1,'K'), comment="""Estimated using an average for rate rule [cCsCJ;CsJ-HH;C]
Euclidian distance = 0
family: 1,2_shiftC"""),
)

reaction(
    label = 'reaction27',
    reactants = ['[O]OC([O])C[C]([O])[O](8116)'],
    products = ['[O]O(16)', '[O][CH]C=C([O])[O](3234)'],
    transitionState = 'TS27',
    kinetics = Arrhenius(A=(1.91635e+11,'s^-1'), n=0.573333, Ea=(139.369,'kJ/mol'), T0=(1,'K'), comment="""Estimated using an average for rate rule [R2OO]
Euclidian distance = 0
Multiplied by reaction path degeneracy 2.0
family: HO2_Elimination_from_PeroxyRadical"""),
)

reaction(
    label = 'reaction28',
    reactants = ['[O]OC([O])C[C]([O])[O](8116)'],
    products = ['[O][C]([O])CC1OOO1(9139)'],
    transitionState = 'TS28',
    kinetics = Arrhenius(A=(1.62e+12,'s^-1'), n=-0.305, Ea=(8.28432,'kJ/mol'), T0=(1,'K'), comment="""Estimated using template [R4_SSS;Y_rad_out;Ypri_rad_out] for rate rule [R4_SSS;O_rad;Opri_rad]
Euclidian distance = 1.41421356237
family: Birad_recombination"""),
)

reaction(
    label = 'reaction29',
    reactants = ['[O]OC([O])C[C]([O])[O](8116)'],
    products = ['[O]OC1CC([O])([O])O1(8118)'],
    transitionState = 'TS29',
    kinetics = Arrhenius(A=(1.62e+12,'s^-1'), n=-0.305, Ea=(8.28432,'kJ/mol'), T0=(1,'K'), comment="""Estimated using template [R4_SSS;Y_rad_out;Ypri_rad_out] for rate rule [R4_SSS;Y_rad_out;Opri_rad]
Euclidian distance = 1.0
family: Birad_recombination"""),
)

reaction(
    label = 'reaction30',
    reactants = ['[O]OC([O])C[C]([O])[O](8116)'],
    products = ['[O]C1CC([O])([O])OO1(9140)'],
    transitionState = 'TS30',
    kinetics = Arrhenius(A=(7.76e+09,'s^-1'), n=0.311, Ea=(7.1128,'kJ/mol'), T0=(1,'K'), comment="""Estimated using template [R5_SSSS;Y_rad_out;Ypri_rad_out] for rate rule [R5_SSSS;O_rad;Ypri_rad_out]
Euclidian distance = 1.0
family: Birad_recombination"""),
)

reaction(
    label = 'reaction31',
    reactants = ['[O]OC([O])C[C]([O])[O](8116)'],
    products = ['[O]OC([O])CC([O])=O(3303)'],
    transitionState = 'TS31',
    kinetics = Arrhenius(A=(1e+10,'s^-1'), n=0, Ea=(0,'kJ/mol'), T0=(1,'K'), comment="""Estimated using an average for rate rule [NOS]
Euclidian distance = 0
family: 1,2-Birad_to_alkene"""),
)

reaction(
    label = 'reaction32',
    reactants = ['O(T)(63)', '[O][C]([O])CC([O])[O](9141)'],
    products = ['[O]OC([O])C[C]([O])[O](8116)'],
    transitionState = 'TS32',
    kinetics = Arrhenius(A=(87.1677,'m^3/(mol*s)'), n=1.88017, Ea=(5.1666,'kJ/mol'), T0=(1,'K'), Tmin=(303.03,'K'), Tmax=(2000,'K'), comment="""From training reaction 3 used for O_rad/NonDe;O_birad
Exact match found for rate rule [O_rad/NonDe;O_birad]
Euclidian distance = 0
Multiplied by reaction path degeneracy 2.0
family: Birad_R_Recombination"""),
)

reaction(
    label = 'reaction33',
    reactants = ['[CH2]C([O])O[O](607)', '[O][C][O](2319)'],
    products = ['[O]OC([O])C[C]([O])[O](8116)'],
    transitionState = 'TS33',
    kinetics = Arrhenius(A=(1.14854e+06,'m^3/(mol*s)'), n=0.575199, Ea=(34.3157,'kJ/mol'), T0=(1,'K'), comment="""Estimated using template [Y_rad;Birad] for rate rule [C_rad/H2/Cs;Birad]
Euclidian distance = 3.0
family: Birad_R_Recombination"""),
)

reaction(
    label = 'reaction34',
    reactants = ['O(T)(63)', '[O]O[CH]C[C]([O])[O](9142)'],
    products = ['[O]OC([O])C[C]([O])[O](8116)'],
    transitionState = 'TS34',
    kinetics = Arrhenius(A=(93609.6,'m^3/(mol*s)'), n=1.13083, Ea=(163.847,'kJ/mol'), T0=(1,'K'), Tmin=(303.03,'K'), Tmax=(2000,'K'), comment="""Estimated using template [Y_rad;O_birad] for rate rule [C_rad/H/CsO;O_birad]
Euclidian distance = 4.0
family: Birad_R_Recombination"""),
)

network(
    label = '2228',
    isomers = [
        '[O]OC([O])C[C]([O])[O](8116)',
    ],
    reactants = [
        ('C=C([O])[O](1172)', '[O]OC=O(5472)'),
    ],
    bathGas = {
        'N2': 0.25,
        'Ne': 0.25,
        'He': 0.25,
        'Ar': 0.25,
    },
)

pressureDependence(
    label = '2228',
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

