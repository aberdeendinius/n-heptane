species(
    label = '[CH2][CH]CC(O)[C]=C(14200)',
    structure = SMILES('[CH2][CH]CC(O)[C]=C'),
    E0 = (400.645,'kJ/mol'),
    modes = [
        HarmonicOscillator(frequencies=([3615,1277.5,1000,2950,3100,1380,975,1025,1650,3025,407.5,1350,352.5,2750,2850,1437.5,1250,1305,750,350,3000,3100,440,815,1455,1000,1685,370,1380,1390,370,380,2900,435,202.925,803.876,1615.5],'cm^-1')),
        HinderedRotor(inertia=(0.153566,'amu*angstrom^2'), symmetry=1, barrier=(3.55403,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.153566,'amu*angstrom^2'), symmetry=1, barrier=(3.55403,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.153566,'amu*angstrom^2'), symmetry=1, barrier=(3.55403,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.153566,'amu*angstrom^2'), symmetry=1, barrier=(3.55403,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.153566,'amu*angstrom^2'), symmetry=1, barrier=(3.55403,'kJ/mol'), semiclassical=False),
    ],
    spinMultiplicity = 4,
    opticalIsomers = 1,
    molecularWeight = (97.1351,'amu'),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[1.16297,0.0662749,-6.09369e-05,3.41313e-08,-8.39586e-12,48285.4,33.0263], Tmin=(100,'K'), Tmax=(945.088,'K')), NASAPolynomial(coeffs=[7.68045,0.0386907,-1.71572e-05,3.24959e-09,-2.2699e-13,47053.5,1.95125], Tmin=(945.088,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(400.645,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(361.68,'J/(mol*K)'), comment="""Thermo group additivity estimation: group(O2s-CsH) + group(Cs-CsCsHH) + group(Cs-(Cds-Cds)CsOsH) + group(Cs-CsCsHH) + group(Cs-CsHHH) + group(Cds-CdsCsH) + group(Cds-CdsHH) + radical(RCCJC) + radical(RCCJ) + radical(Cds_S)"""),
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
    label = 'C=C=CO(12571)',
    structure = SMILES('C=C=CO'),
    E0 = (-26.0646,'kJ/mol'),
    modes = [
        HarmonicOscillator(frequencies=([2950,3100,1380,975,1025,1650,540,610,2055,3615,1277.5,1000,3010,987.5,1337.5,450,1655],'cm^-1')),
        HinderedRotor(inertia=(1.34368,'amu*angstrom^2'), symmetry=1, barrier=(30.8938,'kJ/mol'), semiclassical=False),
    ],
    spinMultiplicity = 1,
    opticalIsomers = 1,
    molecularWeight = (56.0633,'amu'),
    collisionModel = TransportData(shapeIndex=2, epsilon=(3437.21,'J/mol'), sigma=(5.57865,'angstroms'), dipoleMoment=(0,'C*m'), polarizability=(0,'angstroms^3'), rotrelaxcollnum=0, comment="""Epsilon & sigma estimated with Tc=536.88 K, Pc=44.92 bar (from Joback method)"""),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[2.31583,0.0236137,2.05754e-05,-5.73733e-08,2.79863e-11,-3061.58,12.125], Tmin=(100,'K'), Tmax=(901.949,'K')), NASAPolynomial(coeffs=[16.2977,-0.00239911,3.975e-06,-8.57293e-10,5.72973e-14,-7047.88,-62.0029], Tmin=(901.949,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(-26.0646,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(178.761,'J/(mol*K)'), comment="""Thermo group additivity estimation: group(O2s-(Cds-Cd)H) + group(Cds-CdsOsH) + group(Cds-CdsHH) + group(Cdd-CdsCds)"""),
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
    label = '[CH2][CH]CC(O)=C=C(28542)',
    structure = SMILES('[CH2][CH]CC(O)=C=C'),
    E0 = (285.411,'kJ/mol'),
    modes = [
        HarmonicOscillator(frequencies=([3615,1277.5,1000,2950,3100,1380,975,1025,1650,3025,407.5,1350,352.5,350,440,435,1725,540,610,2055,2750,2850,1437.5,1250,1305,750,350,3000,3100,440,815,1455,1000,180,180],'cm^-1')),
        HinderedRotor(inertia=(0.0747035,'amu*angstrom^2'), symmetry=1, barrier=(18.6176,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.0968786,'amu*angstrom^2'), symmetry=1, barrier=(2.22743,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.00338168,'amu*angstrom^2'), symmetry=1, barrier=(18.7041,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.812005,'amu*angstrom^2'), symmetry=1, barrier=(18.6696,'kJ/mol'), semiclassical=False),
    ],
    spinMultiplicity = 3,
    opticalIsomers = 1,
    molecularWeight = (96.1271,'amu'),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[0.605643,0.0702665,-6.62155e-05,3.32622e-08,-6.71257e-12,34452.9,29.9451], Tmin=(100,'K'), Tmax=(1196.3,'K')), NASAPolynomial(coeffs=[14.1205,0.0250767,-9.55214e-06,1.68452e-09,-1.13373e-13,31219.4,-37.6781], Tmin=(1196.3,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(285.411,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(340.893,'J/(mol*K)'), comment="""Thermo group additivity estimation: group(O2s-(Cds-Cd)H) + group(Cs-CsCsHH) + group(Cs-(Cds-Cds)CsHH) + group(Cs-CsHHH) + group(Cds-CdsCsOs) + group(Cds-CdsHH) + group(Cdd-CdsCds) + radical(RCCJC) + radical(RCCJ)"""),
)

species(
    label = 'C=[C]C(O)[CH]C=C(28543)',
    structure = SMILES('C=[C]C(O)[CH]C=C'),
    E0 = (198.879,'kJ/mol'),
    modes = [
        HarmonicOscillator(frequencies=([3615,1277.5,1000,3025,407.5,1350,352.5,2950,3000,3050,3100,1330,1430,900,1050,1000,1050,1600,1700,1380,1390,370,380,2900,435,1685,370,3010,987.5,1337.5,450,1655,434.154,434.329,434.517],'cm^-1')),
        HinderedRotor(inertia=(0.000893809,'amu*angstrom^2'), symmetry=1, barrier=(0.119627,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.191936,'amu*angstrom^2'), symmetry=1, barrier=(25.6876,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.19182,'amu*angstrom^2'), symmetry=1, barrier=(25.6876,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.191976,'amu*angstrom^2'), symmetry=1, barrier=(25.6887,'kJ/mol'), semiclassical=False),
    ],
    spinMultiplicity = 3,
    opticalIsomers = 1,
    molecularWeight = (96.1271,'amu'),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[1.04759,0.0487725,8.83583e-06,-5.00363e-08,2.33403e-11,24040.3,30.1774], Tmin=(100,'K'), Tmax=(978.508,'K')), NASAPolynomial(coeffs=[17.3029,0.0209039,-7.58551e-06,1.44596e-09,-1.07751e-13,19012.1,-57.3297], Tmin=(978.508,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(198.879,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(340.893,'J/(mol*K)'), comment="""Thermo group additivity estimation: group(O2s-CsH) + group(Cs-(Cds-Cds)CsOsH) + group(Cs-(Cds-Cds)CsHH) + group(Cds-CdsCsH) + group(Cds-CdsCsH) + group(Cds-CdsHH) + group(Cds-CdsHH) + radical(Cds_S) + radical(C=CCJC(O)C=C)"""),
)

species(
    label = 'C#CC(O)C[CH][CH2](22868)',
    structure = SMILES('C#CC(O)C[CH][CH2]'),
    E0 = (328.754,'kJ/mol'),
    modes = [
        HarmonicOscillator(frequencies=([750,770,3400,2100,3615,1277.5,1000,3025,407.5,1350,352.5,2750,2850,1437.5,1250,1305,750,350,3000,3100,440,815,1455,1000,2175,525,1380,1390,370,380,2900,435,180,2525.2],'cm^-1')),
        HinderedRotor(inertia=(0.0414107,'amu*angstrom^2'), symmetry=1, barrier=(16.1792,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.0229585,'amu*angstrom^2'), symmetry=1, barrier=(4.18681,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.0169354,'amu*angstrom^2'), symmetry=1, barrier=(76.6546,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.703638,'amu*angstrom^2'), symmetry=1, barrier=(16.178,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(3.33427,'amu*angstrom^2'), symmetry=1, barrier=(76.6615,'kJ/mol'), semiclassical=False),
    ],
    spinMultiplicity = 3,
    opticalIsomers = 1,
    molecularWeight = (96.1271,'amu'),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[0.914362,0.0677339,-6.69348e-05,3.79749e-08,-8.8943e-12,39651.2,31.0509], Tmin=(100,'K'), Tmax=(1024.99,'K')), NASAPolynomial(coeffs=[10.578,0.0300223,-1.17477e-05,2.08113e-09,-1.39794e-13,37670.1,-15.8094], Tmin=(1024.99,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(328.754,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(336.736,'J/(mol*K)'), comment="""Thermo group additivity estimation: group(O2s-CsH) + group(Cs-CsCsHH) + group(Cs-CsCsHH) + group(Cs-(Cds-Cds)CsOsH) + group(Cs-CsHHH) + group(Ct-CtCs) + group(Ct-CtH) + radical(RCCJC) + radical(RCCJ)"""),
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
    label = '[CH2][C]=CO(18753)',
    structure = SMILES('[CH2][C]=CO'),
    E0 = (186.672,'kJ/mol'),
    modes = [
        HarmonicOscillator(frequencies=([3000,3100,440,815,1455,1000,1685,370,3615,1277.5,1000,3010,987.5,1337.5,450,1655],'cm^-1')),
        HinderedRotor(inertia=(1.23523,'amu*angstrom^2'), symmetry=1, barrier=(28.4004,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(1.2351,'amu*angstrom^2'), symmetry=1, barrier=(28.3973,'kJ/mol'), semiclassical=False),
    ],
    spinMultiplicity = 3,
    opticalIsomers = 1,
    molecularWeight = (56.0633,'amu'),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[2.24497,0.0260528,1.3484e-05,-5.00525e-08,2.54383e-11,22526.5,14.0801], Tmin=(100,'K'), Tmax=(898.827,'K')), NASAPolynomial(coeffs=[16.2027,-0.00210248,3.79693e-06,-8.3211e-10,5.63273e-14,18645.6,-59.4], Tmin=(898.827,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(186.672,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(174.604,'J/(mol*K)'), comment="""Thermo group additivity estimation: group(O2s-(Cds-Cd)H) + group(Cs-(Cds-Cds)HHH) + group(Cds-CdsCsH) + group(Cds-CdsOsH) + radical(Cds_S) + radical(Allyl_P)"""),
)

species(
    label = 'OH(D)(132)',
    structure = SMILES('[OH]'),
    E0 = (28.3945,'kJ/mol'),
    modes = [
        HarmonicOscillator(frequencies=([3668.68],'cm^-1')),
    ],
    spinMultiplicity = 2,
    opticalIsomers = 1,
    molecularWeight = (17.0073,'amu'),
    collisionModel = TransportData(shapeIndex=1, epsilon=(665.16,'J/mol'), sigma=(2.75,'angstroms'), dipoleMoment=(0,'C*m'), polarizability=(0,'angstroms^3'), rotrelaxcollnum=0.0, comment="""GRI-Mech"""),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[3.51457,2.92814e-05,-5.32177e-07,1.01951e-09,-3.85951e-13,3414.25,2.10435], Tmin=(100,'K'), Tmax=(1145.75,'K')), NASAPolynomial(coeffs=[3.07194,0.000604011,-1.39759e-08,-2.13452e-11,2.4807e-15,3579.39,4.57799], Tmin=(1145.75,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(28.3945,'kJ/mol'), Cp0=(29.1007,'J/(mol*K)'), CpInf=(37.4151,'J/(mol*K)'), label="""OH(D)""", comment="""Thermo library: primaryThermoLibrary"""),
)

species(
    label = '[CH2][CH]CC=C=C(17733)',
    structure = SMILES('[CH2][CH]CC=C=C'),
    E0 = (499.969,'kJ/mol'),
    modes = [
        HarmonicOscillator(frequencies=([2950,3100,1380,975,1025,1650,3025,407.5,1350,352.5,540,610,2055,2750,2850,1437.5,1250,1305,750,350,3000,3100,440,815,1455,1000,3010,987.5,1337.5,450,1655,180,757.009],'cm^-1')),
        HinderedRotor(inertia=(0.106262,'amu*angstrom^2'), symmetry=1, barrier=(2.44316,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.106185,'amu*angstrom^2'), symmetry=1, barrier=(2.4414,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.106216,'amu*angstrom^2'), symmetry=1, barrier=(2.44212,'kJ/mol'), semiclassical=False),
    ],
    spinMultiplicity = 3,
    opticalIsomers = 1,
    molecularWeight = (80.1277,'amu'),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[1.99087,0.047132,-2.8655e-05,9.31653e-09,-1.34472e-12,60201.8,24.9057], Tmin=(100,'K'), Tmax=(1434.16,'K')), NASAPolynomial(coeffs=[6.91729,0.0333917,-1.4284e-05,2.63616e-09,-1.8021e-13,58788.7,-0.637827], Tmin=(1434.16,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(499.969,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(320.107,'J/(mol*K)'), comment="""Thermo group additivity estimation: group(Cs-CsCsHH) + group(Cs-(Cds-Cds)CsHH) + group(Cs-CsHHH) + group(Cds-CdsCsH) + group(Cds-CdsHH) + group(Cdd-CdsCds) + radical(RCCJ) + radical(RCCJC)"""),
)

species(
    label = '[CH2][CH]CC(O)=C[CH2](14197)',
    structure = SMILES('[CH2][CH]CC(O)=C[CH2]'),
    E0 = (260.306,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
    molecularWeight = (97.1351,'amu'),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[0.486036,0.0690806,-5.20013e-05,1.59787e-08,-3.67205e-13,31441.2,31.6008], Tmin=(100,'K'), Tmax=(1027.75,'K')), NASAPolynomial(coeffs=[15.1728,0.0259344,-9.48405e-06,1.66754e-09,-1.13548e-13,27682.1,-43.2574], Tmin=(1027.75,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(260.306,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(361.68,'J/(mol*K)'), comment="""Thermo group additivity estimation: group(O2s-(Cds-Cd)H) + group(Cs-CsCsHH) + group(Cs-(Cds-Cds)CsHH) + group(Cs-CsHHH) + group(Cs-(Cds-Cds)HHH) + group(Cds-CdsCsOs) + group(Cds-CdsCsH) + radical(Allyl_P) + radical(RCCJC) + radical(RCCJ)"""),
)

species(
    label = '[CH2]C[CH]C(O)[C]=C(28544)',
    structure = SMILES('[CH2]C[CH]C(O)[C]=C'),
    E0 = (406.101,'kJ/mol'),
    modes = [
        HarmonicOscillator(frequencies=([3615,1277.5,1000,2950,3100,1380,975,1025,1650,3025,407.5,1350,352.5,2750,2850,1437.5,1250,1305,750,350,3000,3100,440,815,1455,1000,1685,370,1380,1390,370,380,2900,435,200,800,1600],'cm^-1')),
        HinderedRotor(inertia=(0.156089,'amu*angstrom^2'), symmetry=1, barrier=(3.5888,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.156089,'amu*angstrom^2'), symmetry=1, barrier=(3.5888,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.156089,'amu*angstrom^2'), symmetry=1, barrier=(3.5888,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.156089,'amu*angstrom^2'), symmetry=1, barrier=(3.5888,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.156089,'amu*angstrom^2'), symmetry=1, barrier=(3.5888,'kJ/mol'), semiclassical=False),
    ],
    spinMultiplicity = 4,
    opticalIsomers = 1,
    molecularWeight = (97.1351,'amu'),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[0.697355,0.0663869,-5.24802e-05,2.15784e-08,-3.60183e-12,48966.7,34.4496], Tmin=(100,'K'), Tmax=(1413.76,'K')), NASAPolynomial(coeffs=[14.5546,0.0271798,-1.08812e-05,1.96205e-09,-1.32984e-13,45048.5,-37.2017], Tmin=(1413.76,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(406.101,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(361.68,'J/(mol*K)'), comment="""Thermo group additivity estimation: group(O2s-CsH) + group(Cs-CsCsHH) + group(Cs-(Cds-Cds)CsOsH) + group(Cs-CsCsHH) + group(Cs-CsHHH) + group(Cds-CdsCsH) + group(Cds-CdsHH) + radical(CCJCO) + radical(RCCJ) + radical(Cds_S)"""),
)

species(
    label = '[CH]=CC(O)C[CH][CH2](14203)',
    structure = SMILES('[CH]=CC(O)C[CH][CH2]'),
    E0 = (409.9,'kJ/mol'),
    modes = [
        HarmonicOscillator(frequencies=([3120,650,792.5,1650,3615,1277.5,1000,3025,407.5,1350,352.5,2750,2850,1437.5,1250,1305,750,350,3000,3100,440,815,1455,1000,3010,987.5,1337.5,450,1655,1380,1390,370,380,2900,435,1906.25,2423.01],'cm^-1')),
        HinderedRotor(inertia=(0.086439,'amu*angstrom^2'), symmetry=1, barrier=(6.39769,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.260917,'amu*angstrom^2'), symmetry=1, barrier=(19.3002,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.0864225,'amu*angstrom^2'), symmetry=1, barrier=(6.39756,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.0016165,'amu*angstrom^2'), symmetry=1, barrier=(0.119627,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.260867,'amu*angstrom^2'), symmetry=1, barrier=(19.2995,'kJ/mol'), semiclassical=False),
    ],
    spinMultiplicity = 4,
    opticalIsomers = 1,
    molecularWeight = (97.1351,'amu'),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[1.0747,0.0658921,-5.57874e-05,2.67338e-08,-5.46751e-12,49403.6,33.2627], Tmin=(100,'K'), Tmax=(1133.58,'K')), NASAPolynomial(coeffs=[9.77838,0.0351802,-1.51484e-05,2.8339e-09,-1.9666e-13,47430.3,-9.81881], Tmin=(1133.58,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(409.9,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(361.68,'J/(mol*K)'), comment="""Thermo group additivity estimation: group(O2s-CsH) + group(Cs-CsCsHH) + group(Cs-(Cds-Cds)CsOsH) + group(Cs-CsCsHH) + group(Cs-CsHHH) + group(Cds-CdsCsH) + group(Cds-CdsHH) + radical(RCCJC) + radical(RCCJ) + radical(Cds_P)"""),
)

species(
    label = '[CH2][CH]CC([O])C=C(12801)',
    structure = SMILES('[CH2][CH]CC([O])C=C'),
    E0 = (393.164,'kJ/mol'),
    modes = [
        HarmonicOscillator(frequencies=([2950,3100,1380,975,1025,1650,3025,407.5,1350,352.5,2750,2850,1437.5,1250,1305,750,350,3000,3100,440,815,1455,1000,3010,987.5,1337.5,450,1655,1380,1390,370,380,2900,435,335.453,1351.14,1360.57,1364.23],'cm^-1')),
        HinderedRotor(inertia=(0.0915869,'amu*angstrom^2'), symmetry=1, barrier=(7.3102,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.0888427,'amu*angstrom^2'), symmetry=1, barrier=(7.30888,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.00210686,'amu*angstrom^2'), symmetry=1, barrier=(0.171365,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.232374,'amu*angstrom^2'), symmetry=1, barrier=(18.5645,'kJ/mol'), semiclassical=False),
    ],
    spinMultiplicity = 4,
    opticalIsomers = 1,
    molecularWeight = (97.1351,'amu'),
    collisionModel = TransportData(shapeIndex=2, epsilon=(3797.28,'J/mol'), sigma=(6.5936,'angstroms'), dipoleMoment=(0,'C*m'), polarizability=(0,'angstroms^3'), rotrelaxcollnum=0, comment="""Epsilon & sigma estimated with Tc=593.13 K, Pc=30.06 bar (from Joback method)"""),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[1.11089,0.0607652,-4.15362e-05,1.4774e-08,-2.19158e-12,47392.6,32.5333], Tmin=(100,'K'), Tmax=(1523.24,'K')), NASAPolynomial(coeffs=[12.028,0.032097,-1.33054e-05,2.41832e-09,-1.63715e-13,44066.8,-24.7298], Tmin=(1523.24,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(393.164,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(365.837,'J/(mol*K)'), comment="""Thermo group additivity estimation: group(O2s-CsH) + group(Cs-CsCsHH) + group(Cs-(Cds-Cds)CsOsH) + group(Cs-CsCsHH) + group(Cs-CsHHH) + group(Cds-CdsCsH) + group(Cds-CdsHH) + radical(CC(C)OJ) + radical(RCCJ) + radical(RCCJC)"""),
)

species(
    label = '[CH2][CH][CH]C(O)C=C(14199)',
    structure = SMILES('[CH2][CH][CH]C(O)C=C'),
    E0 = (362.706,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
    molecularWeight = (97.1351,'amu'),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[1.02105,0.060939,-4.35701e-05,1.64184e-08,-2.56116e-12,43734.2,35.5475], Tmin=(100,'K'), Tmax=(1476.65,'K')), NASAPolynomial(coeffs=[12.5126,0.0298103,-1.19493e-05,2.14244e-09,-1.44213e-13,40340.4,-24.3717], Tmin=(1476.65,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(362.706,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(361.68,'J/(mol*K)'), comment="""Thermo group additivity estimation: group(O2s-CsH) + group(Cs-CsCsHH) + group(Cs-(Cds-Cds)CsOsH) + group(Cs-CsCsHH) + group(Cs-CsHHH) + group(Cds-CdsCsH) + group(Cds-CdsHH) + radical(CCJCO) + radical(RCCJC) + radical(RCCJ)"""),
)

species(
    label = '[CH2][C]=C(O)CC[CH2](28545)',
    structure = SMILES('[CH2][C]=C(O)CC[CH2]'),
    E0 = (303.701,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
    molecularWeight = (97.1351,'amu'),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[0.233242,0.0737237,-5.82199e-05,1.77929e-08,-3.47969e-14,36670.4,30.2462], Tmin=(100,'K'), Tmax=(1009.6,'K')), NASAPolynomial(coeffs=[17.2198,0.0233186,-8.43287e-06,1.49247e-09,-1.02838e-13,32379.4,-56.1309], Tmin=(1009.6,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(303.701,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(361.68,'J/(mol*K)'), comment="""Thermo group additivity estimation: group(O2s-(Cds-Cd)H) + group(Cs-CsCsHH) + group(Cs-(Cds-Cds)CsHH) + group(Cs-CsHHH) + group(Cs-(Cds-Cds)HHH) + group(Cds-CdsCsOs) + group(Cds-CdsCsH) + radical(Cds_S) + radical(Allyl_P) + radical(RCCJ)"""),
)

species(
    label = 'C=[C]C(O)[CH][CH]C(27904)',
    structure = SMILES('C=[C]C(O)[CH][CH]C'),
    E0 = (395.301,'kJ/mol'),
    modes = [
        HarmonicOscillator(frequencies=([3615,1277.5,1000,3000,3050,390,425,1340,1360,335,370,2950,3100,1380,975,1025,1650,1380,1390,370,380,2900,435,2750,2800,2850,1350,1500,750,1050,1375,1000,1685,370,232.809,874.638,1862.47],'cm^-1')),
        HinderedRotor(inertia=(0.106337,'amu*angstrom^2'), symmetry=1, barrier=(3.31283,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.106337,'amu*angstrom^2'), symmetry=1, barrier=(3.31283,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.106337,'amu*angstrom^2'), symmetry=1, barrier=(3.31283,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.106337,'amu*angstrom^2'), symmetry=1, barrier=(3.31283,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.106337,'amu*angstrom^2'), symmetry=1, barrier=(3.31283,'kJ/mol'), semiclassical=False),
    ],
    spinMultiplicity = 4,
    opticalIsomers = 1,
    molecularWeight = (97.1351,'amu'),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[1.22725,0.0616182,-4.71496e-05,2.00992e-08,-3.67605e-12,47642.9,33.491], Tmin=(100,'K'), Tmax=(1245.38,'K')), NASAPolynomial(coeffs=[9.61903,0.0346648,-1.46855e-05,2.72078e-09,-1.87461e-13,45552.7,-8.83598], Tmin=(1245.38,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(395.301,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(361.68,'J/(mol*K)'), comment="""Thermo group additivity estimation: group(O2s-CsH) + group(Cs-CsCsHH) + group(Cs-(Cds-Cds)CsOsH) + group(Cs-CsCsHH) + group(Cs-CsHHH) + group(Cds-CdsCsH) + group(Cds-CdsHH) + radical(CCJCO) + radical(RCCJC) + radical(Cds_S)"""),
)

species(
    label = '[CH2]CCC([O])[C]=C(14202)',
    structure = SMILES('[CH2]CCC([O])[C]=C'),
    E0 = (436.56,'kJ/mol'),
    modes = [
        HarmonicOscillator(frequencies=([1685,370,2750,2783.33,2816.67,2850,1425,1450,1225,1275,1270,1340,700,800,300,400,2950,3100,1380,975,1025,1650,1380,1390,370,380,2900,435,3000,3100,440,815,1455,1000,180,180,999.421,4000],'cm^-1')),
        HinderedRotor(inertia=(0.838523,'amu*angstrom^2'), symmetry=1, barrier=(19.2793,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.0349809,'amu*angstrom^2'), symmetry=1, barrier=(2.99919,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.0271975,'amu*angstrom^2'), symmetry=1, barrier=(19.2788,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.838507,'amu*angstrom^2'), symmetry=1, barrier=(19.2789,'kJ/mol'), semiclassical=False),
    ],
    spinMultiplicity = 4,
    opticalIsomers = 1,
    molecularWeight = (97.1351,'amu'),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[0.784842,0.0662291,-5.04654e-05,1.99314e-08,-3.22619e-12,52625.2,31.4449], Tmin=(100,'K'), Tmax=(1437.76,'K')), NASAPolynomial(coeffs=[13.9091,0.0297163,-1.23722e-05,2.2683e-09,-1.54915e-13,48851.3,-36.6371], Tmin=(1437.76,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(436.56,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(365.837,'J/(mol*K)'), comment="""Thermo group additivity estimation: group(O2s-CsH) + group(Cs-CsCsHH) + group(Cs-(Cds-Cds)CsOsH) + group(Cs-CsCsHH) + group(Cs-CsHHH) + group(Cds-CdsCsH) + group(Cds-CdsHH) + radical(RCCJ) + radical(CC(C)OJ) + radical(Cds_S)"""),
)

species(
    label = '[CH2][C]=C(O)C[CH]C(27905)',
    structure = SMILES('[CH2][C]=C(O)C[CH]C'),
    E0 = (292.901,'kJ/mol'),
    modes = [
        HarmonicOscillator(frequencies=([3615,1277.5,1000,3025,407.5,1350,352.5,350,440,435,1725,2750,2850,1437.5,1250,1305,750,350,3000,3100,440,815,1455,1000,2750,2800,2850,1350,1500,750,1050,1375,1000,1685,370,180,180],'cm^-1')),
        HinderedRotor(inertia=(0.000520282,'amu*angstrom^2'), symmetry=1, barrier=(5.73083,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.0344301,'amu*angstrom^2'), symmetry=1, barrier=(17.0028,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.739537,'amu*angstrom^2'), symmetry=1, barrier=(17.0034,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.034436,'amu*angstrom^2'), symmetry=1, barrier=(17.0045,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(3.81122,'amu*angstrom^2'), symmetry=1, barrier=(87.6274,'kJ/mol'), semiclassical=False),
    ],
    spinMultiplicity = 4,
    opticalIsomers = 1,
    molecularWeight = (97.1351,'amu'),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[0.438077,0.0726667,-6.52427e-05,3.1349e-08,-6.0738e-12,35360.8,30.461], Tmin=(100,'K'), Tmax=(1241.51,'K')), NASAPolynomial(coeffs=[14.585,0.0270863,-1.01713e-05,1.77634e-09,-1.18734e-13,31848.2,-40.85], Tmin=(1241.51,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(292.901,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(361.68,'J/(mol*K)'), comment="""Thermo group additivity estimation: group(O2s-(Cds-Cd)H) + group(Cs-CsCsHH) + group(Cs-(Cds-Cds)CsHH) + group(Cs-CsHHH) + group(Cs-(Cds-Cds)HHH) + group(Cds-CdsCsOs) + group(Cds-CdsCsH) + radical(RCCJC) + radical(Allyl_P) + radical(Cds_S)"""),
)

species(
    label = '[CH]=[C]C(O)CC[CH2](28546)',
    structure = SMILES('[CH]=[C]C(O)CC[CH2]'),
    E0 = (453.295,'kJ/mol'),
    modes = [
        HarmonicOscillator(frequencies=([3120,650,792.5,1650,3615,1277.5,1000,2750,2783.33,2816.67,2850,1425,1450,1225,1275,1270,1340,700,800,300,400,3000,3100,440,815,1455,1000,1685,370,1380,1390,370,380,2900,435,276.125,276.151],'cm^-1')),
        HinderedRotor(inertia=(0.0022066,'amu*angstrom^2'), symmetry=1, barrier=(0.119627,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.00220933,'amu*angstrom^2'), symmetry=1, barrier=(0.119627,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.00220981,'amu*angstrom^2'), symmetry=1, barrier=(0.119627,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.208439,'amu*angstrom^2'), symmetry=1, barrier=(11.2748,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.208305,'amu*angstrom^2'), symmetry=1, barrier=(11.2747,'kJ/mol'), semiclassical=False),
    ],
    spinMultiplicity = 4,
    opticalIsomers = 1,
    molecularWeight = (97.1351,'amu'),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[0.766314,0.0712011,-6.43616e-05,3.16008e-08,-6.42207e-12,54635.2,32.1065], Tmin=(100,'K'), Tmax=(1162.07,'K')), NASAPolynomial(coeffs=[12.1496,0.0320184,-1.37845e-05,2.58531e-09,-1.79874e-13,51989.5,-24.521], Tmin=(1162.07,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(453.295,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(361.68,'J/(mol*K)'), comment="""Thermo group additivity estimation: group(O2s-CsH) + group(Cs-CsCsHH) + group(Cs-(Cds-Cds)CsOsH) + group(Cs-CsCsHH) + group(Cs-CsHHH) + group(Cds-CdsCsH) + group(Cds-CdsHH) + radical(RCCJ) + radical(Cds_P) + radical(Cds_S)"""),
)

species(
    label = 'C=[C]C([O])C[CH]C(13187)',
    structure = SMILES('C=[C]C([O])C[CH]C'),
    E0 = (425.76,'kJ/mol'),
    modes = [
        HarmonicOscillator(frequencies=([2950,3100,1380,975,1025,1650,3025,407.5,1350,352.5,2750,2850,1437.5,1250,1305,750,350,1380,1390,370,380,2900,435,2750,2800,2850,1350,1500,750,1050,1375,1000,1685,370,445.2,445.2,445.2,3551.95],'cm^-1')),
        HinderedRotor(inertia=(0.0601108,'amu*angstrom^2'), symmetry=1, barrier=(8.45448,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.0601105,'amu*angstrom^2'), symmetry=1, barrier=(8.45449,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.0601105,'amu*angstrom^2'), symmetry=1, barrier=(8.45448,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.0280813,'amu*angstrom^2'), symmetry=1, barrier=(27.1168,'kJ/mol'), semiclassical=False),
    ],
    spinMultiplicity = 4,
    opticalIsomers = 1,
    molecularWeight = (97.1351,'amu'),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[1.30749,0.0615208,-4.5259e-05,1.85347e-08,-3.31578e-12,51301.9,30.5144], Tmin=(100,'K'), Tmax=(1246.25,'K')), NASAPolynomial(coeffs=[8.79305,0.0374947,-1.63405e-05,3.06496e-09,-2.12489e-13,49436.1,-7.24688], Tmin=(1246.25,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(425.76,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(365.837,'J/(mol*K)'), comment="""Thermo group additivity estimation: group(O2s-CsH) + group(Cs-CsCsHH) + group(Cs-(Cds-Cds)CsOsH) + group(Cs-CsCsHH) + group(Cs-CsHHH) + group(Cds-CdsCsH) + group(Cds-CdsHH) + radical(RCCJC) + radical(Cds_S) + radical(CC(C)OJ)"""),
)

species(
    label = '[CH]=[C]C(O)C[CH]C(27906)',
    structure = SMILES('[CH]=[C]C(O)C[CH]C'),
    E0 = (442.495,'kJ/mol'),
    modes = [
        HarmonicOscillator(frequencies=([3120,650,792.5,1650,3615,1277.5,1000,3025,407.5,1350,352.5,2750,2850,1437.5,1250,1305,750,350,1380,1390,370,380,2900,435,2750,2800,2850,1350,1500,750,1050,1375,1000,1685,370,180,2400.99],'cm^-1')),
        HinderedRotor(inertia=(0.118148,'amu*angstrom^2'), symmetry=1, barrier=(2.71646,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.452344,'amu*angstrom^2'), symmetry=1, barrier=(10.4003,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.0363792,'amu*angstrom^2'), symmetry=1, barrier=(10.4003,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.036379,'amu*angstrom^2'), symmetry=1, barrier=(10.4003,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.0363788,'amu*angstrom^2'), symmetry=1, barrier=(10.4003,'kJ/mol'), semiclassical=False),
    ],
    spinMultiplicity = 4,
    opticalIsomers = 1,
    molecularWeight = (97.1351,'amu'),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[1.022,0.0694113,-6.83306e-05,4.07754e-08,-1.05157e-11,53323.8,32.1493], Tmin=(100,'K'), Tmax=(912.955,'K')), NASAPolynomial(coeffs=[8.17293,0.0380803,-1.68531e-05,3.18485e-09,-2.22053e-13,52018.1,-1.69858], Tmin=(912.955,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(442.495,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(361.68,'J/(mol*K)'), comment="""Thermo group additivity estimation: group(O2s-CsH) + group(Cs-CsCsHH) + group(Cs-(Cds-Cds)CsOsH) + group(Cs-CsCsHH) + group(Cs-CsHHH) + group(Cds-CdsCsH) + group(Cds-CdsHH) + radical(RCCJC) + radical(Cds_P) + radical(Cds_S)"""),
)

species(
    label = '[CH2][C]=CC[CH][CH2](17739)',
    structure = SMILES('[CH2][C]=CC[CH][CH2]'),
    E0 = (712.706,'kJ/mol'),
    modes = [
        HarmonicOscillator(frequencies=([3000,3033.33,3066.67,3100,415,465,780,850,1435,1475,900,1100,2750,2850,1437.5,1250,1305,750,350,3010,987.5,1337.5,450,1655,1685,370,3025,407.5,1350,352.5,261.72,1585.42],'cm^-1')),
        HinderedRotor(inertia=(0.0024611,'amu*angstrom^2'), symmetry=1, barrier=(0.11981,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.0961317,'amu*angstrom^2'), symmetry=1, barrier=(4.66993,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.0356867,'amu*angstrom^2'), symmetry=1, barrier=(63.6449,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.096257,'amu*angstrom^2'), symmetry=1, barrier=(4.67159,'kJ/mol'), semiclassical=False),
    ],
    spinMultiplicity = 5,
    opticalIsomers = 1,
    molecularWeight = (80.1277,'amu'),
    collisionModel = TransportData(shapeIndex=2, epsilon=(3287.25,'J/mol'), sigma=(5.9252,'angstroms'), dipoleMoment=(0,'C*m'), polarizability=(0,'angstroms^3'), rotrelaxcollnum=0, comment="""Epsilon & sigma estimated with Tc=513.46 K, Pc=35.86 bar (from Joback method)"""),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[1.89726,0.0496265,-3.51675e-05,1.51569e-08,-3.07261e-12,85791.6,26.9596], Tmin=(100,'K'), Tmax=(1064.4,'K')), NASAPolynomial(coeffs=[5.57814,0.0357938,-1.56738e-05,2.94733e-09,-2.04892e-13,85008,8.97174], Tmin=(1064.4,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(712.706,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(315.95,'J/(mol*K)'), comment="""Thermo group additivity estimation: group(Cs-CsCsHH) + group(Cs-(Cds-Cds)CsHH) + group(Cs-CsHHH) + group(Cs-(Cds-Cds)HHH) + group(Cds-CdsCsH) + group(Cds-CdsCsH) + radical(Allyl_P) + radical(RCCJ) + radical(RCCJC) + radical(Cds_S)"""),
)

species(
    label = '[CH2][CH]CC([O])[C]=C(14207)',
    structure = SMILES('[CH2][CH]CC([O])[C]=C'),
    E0 = (631.006,'kJ/mol'),
    modes = [
        HarmonicOscillator(frequencies=([2950,3100,1380,975,1025,1650,3025,407.5,1350,352.5,2750,2850,1437.5,1250,1305,750,350,3000,3100,440,815,1455,1000,1685,370,1380,1390,370,380,2900,435,180,853.153,2641.76,2643.64],'cm^-1')),
        HinderedRotor(inertia=(0.228319,'amu*angstrom^2'), symmetry=1, barrier=(5.24951,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.578186,'amu*angstrom^2'), symmetry=1, barrier=(13.2936,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.064128,'amu*angstrom^2'), symmetry=1, barrier=(34.3497,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.0245751,'amu*angstrom^2'), symmetry=1, barrier=(13.8631,'kJ/mol'), semiclassical=False),
    ],
    spinMultiplicity = 5,
    opticalIsomers = 1,
    molecularWeight = (96.1271,'amu'),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[1.35025,0.0620229,-5.42584e-05,2.82452e-08,-6.47739e-12,75984.8,32.1602], Tmin=(100,'K'), Tmax=(1003,'K')), NASAPolynomial(coeffs=[7.68232,0.0367699,-1.64915e-05,3.14212e-09,-2.20275e-13,74714.6,1.59277], Tmin=(1003,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(631.006,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(340.893,'J/(mol*K)'), comment="""Thermo group additivity estimation: group(O2s-CsH) + group(Cs-CsCsHH) + group(Cs-(Cds-Cds)CsOsH) + group(Cs-CsCsHH) + group(Cs-CsHHH) + group(Cds-CdsCsH) + group(Cds-CdsHH) + radical(CC(C)OJ) + radical(Cds_S) + radical(RCCJC) + radical(RCCJ)"""),
)

species(
    label = '[CH2][C]=C(O)C[CH][CH2](28547)',
    structure = SMILES('[CH2][C]=C(O)C[CH][CH2]'),
    E0 = (498.147,'kJ/mol'),
    modes = [
        HarmonicOscillator(frequencies=([3615,1277.5,1000,3025,407.5,1350,352.5,350,440,435,1725,2750,2850,1437.5,1250,1305,750,350,3000,3033.33,3066.67,3100,415,465,780,850,1435,1475,900,1100,1685,370,180,1185.21],'cm^-1')),
        HinderedRotor(inertia=(0.00250639,'amu*angstrom^2'), symmetry=1, barrier=(2.49675,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.728804,'amu*angstrom^2'), symmetry=1, barrier=(16.7566,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.284911,'amu*angstrom^2'), symmetry=1, barrier=(16.7564,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.729278,'amu*angstrom^2'), symmetry=1, barrier=(16.7675,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.0113687,'amu*angstrom^2'), symmetry=1, barrier=(81.5924,'kJ/mol'), semiclassical=False),
    ],
    spinMultiplicity = 5,
    opticalIsomers = 1,
    molecularWeight = (96.1271,'amu'),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[0.596356,0.0719821,-7.08089e-05,3.74244e-08,-7.96091e-12,60038.3,31.6794], Tmin=(100,'K'), Tmax=(1136.35,'K')), NASAPolynomial(coeffs=[13.6689,0.025966,-1.00665e-05,1.78822e-09,-1.20797e-13,57067.3,-33.0588], Tmin=(1136.35,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(498.147,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(336.736,'J/(mol*K)'), comment="""Thermo group additivity estimation: group(O2s-(Cds-Cd)H) + group(Cs-CsCsHH) + group(Cs-(Cds-Cds)CsHH) + group(Cs-CsHHH) + group(Cs-(Cds-Cds)HHH) + group(Cds-CdsCsOs) + group(Cds-CdsCsH) + radical(Allyl_P) + radical(Cds_S) + radical(RCCJ) + radical(RCCJC)"""),
)

species(
    label = '[CH2][CH][CH]C(O)[C]=C(28548)',
    structure = SMILES('[CH2][CH][CH]C(O)[C]=C'),
    E0 = (600.548,'kJ/mol'),
    modes = [
        HarmonicOscillator(frequencies=([3615,1277.5,1000,3000,3050,390,425,1340,1360,335,370,2950,3100,1380,975,1025,1650,3000,3100,440,815,1455,1000,1685,370,1380,1390,370,380,2900,435,228.73,649.77,3336.81],'cm^-1')),
        HinderedRotor(inertia=(0.0430468,'amu*angstrom^2'), symmetry=1, barrier=(1.48923,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.0430468,'amu*angstrom^2'), symmetry=1, barrier=(1.48923,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.0430468,'amu*angstrom^2'), symmetry=1, barrier=(1.48923,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.0430468,'amu*angstrom^2'), symmetry=1, barrier=(1.48923,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.0430468,'amu*angstrom^2'), symmetry=1, barrier=(1.48923,'kJ/mol'), semiclassical=False),
    ],
    spinMultiplicity = 5,
    opticalIsomers = 1,
    molecularWeight = (96.1271,'amu'),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[1.31053,0.0616994,-5.49123e-05,2.84758e-08,-6.3602e-12,72323.9,34.9873], Tmin=(100,'K'), Tmax=(1041.88,'K')), NASAPolynomial(coeffs=[8.57705,0.0338015,-1.47472e-05,2.77519e-09,-1.93255e-13,70809.7,-0.367525], Tmin=(1041.88,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(600.548,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(336.736,'J/(mol*K)'), comment="""Thermo group additivity estimation: group(O2s-CsH) + group(Cs-CsCsHH) + group(Cs-(Cds-Cds)CsOsH) + group(Cs-CsCsHH) + group(Cs-CsHHH) + group(Cds-CdsCsH) + group(Cds-CdsHH) + radical(RCCJ) + radical(CCJCO) + radical(Cds_S) + radical(RCCJC)"""),
)

species(
    label = '[CH]=[C]C(O)C[CH][CH2](28549)',
    structure = SMILES('[CH]=[C]C(O)C[CH][CH2]'),
    E0 = (647.742,'kJ/mol'),
    modes = [
        HarmonicOscillator(frequencies=([3120,650,792.5,1650,3615,1277.5,1000,3025,407.5,1350,352.5,2750,2850,1437.5,1250,1305,750,350,3000,3100,440,815,1455,1000,1685,370,1380,1390,370,380,2900,435,414.625,1685.67],'cm^-1')),
        HinderedRotor(inertia=(0.000980609,'amu*angstrom^2'), symmetry=1, barrier=(0.119627,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.0537677,'amu*angstrom^2'), symmetry=1, barrier=(6.55887,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.0537724,'amu*angstrom^2'), symmetry=1, barrier=(6.55842,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.0537841,'amu*angstrom^2'), symmetry=1, barrier=(6.5588,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.0537569,'amu*angstrom^2'), symmetry=1, barrier=(6.55915,'kJ/mol'), semiclassical=False),
    ],
    spinMultiplicity = 5,
    opticalIsomers = 1,
    molecularWeight = (96.1271,'amu'),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[0.955446,0.0714151,-8.34644e-05,5.95166e-08,-1.79456e-11,78010.9,34.1715], Tmin=(100,'K'), Tmax=(796.419,'K')), NASAPolynomial(coeffs=[8.09096,0.0355769,-1.59655e-05,3.01438e-09,-2.09155e-13,76874.3,1.37104], Tmin=(796.419,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(647.742,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(336.736,'J/(mol*K)'), comment="""Thermo group additivity estimation: group(O2s-CsH) + group(Cs-CsCsHH) + group(Cs-(Cds-Cds)CsOsH) + group(Cs-CsCsHH) + group(Cs-CsHHH) + group(Cds-CdsCsH) + group(Cds-CdsHH) + radical(RCCJ) + radical(RCCJC) + radical(Cds_P) + radical(Cds_S)"""),
)

species(
    label = '[CH2]CCC(O)=C=C(28550)',
    structure = SMILES('[CH2]CCC(O)=C=C'),
    E0 = (90.9645,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
    molecularWeight = (97.1351,'amu'),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[0.316084,0.0711408,-5.06146e-05,9.79145e-09,2.80998e-12,11081.9,28.2482], Tmin=(100,'K'), Tmax=(998.252,'K')), NASAPolynomial(coeffs=[17.2624,0.0231105,-8.30561e-06,1.47924e-09,-1.02856e-13,6708.32,-58.4377], Tmin=(998.252,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(90.9645,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(365.837,'J/(mol*K)'), comment="""Thermo group additivity estimation: group(O2s-(Cds-Cd)H) + group(Cs-CsCsHH) + group(Cs-(Cds-Cds)CsHH) + group(Cs-CsHHH) + group(Cds-CdsCsOs) + group(Cds-CdsHH) + group(Cdd-CdsCds) + radical(RCCJ)"""),
)

species(
    label = 'C=C[CH]C(O)C=C(14211)',
    structure = SMILES('C=C[CH]C(O)C=C'),
    E0 = (-38.9632,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
    molecularWeight = (97.1351,'amu'),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[1.08192,0.044158,3.3673e-05,-7.93709e-08,3.43363e-11,-4563.09,29.5804], Tmin=(100,'K'), Tmax=(967.453,'K')), NASAPolynomial(coeffs=[18.4602,0.0214577,-7.3386e-06,1.40416e-09,-1.07018e-13,-10225.8,-65.5726], Tmin=(967.453,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(-38.9632,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(365.837,'J/(mol*K)'), comment="""Thermo group additivity estimation: group(O2s-CsH) + group(Cs-(Cds-Cds)CsOsH) + group(Cs-(Cds-Cds)CsHH) + group(Cds-CdsCsH) + group(Cds-CdsCsH) + group(Cds-CdsHH) + group(Cds-CdsHH) + radical(C=CCJC(O)C=C)"""),
)

species(
    label = '[CH2][CH]CC([CH2])=CO(14139)',
    structure = SMILES('[CH2][CH]CC([CH2])=CO'),
    E0 = (263.042,'kJ/mol'),
    modes = [
        HarmonicOscillator(frequencies=([3615,1277.5,1000,3025,407.5,1350,352.5,350,440,435,1725,2750,2850,1437.5,1250,1305,750,350,3000,3033.33,3066.67,3100,415,465,780,850,1435,1475,900,1100,3010,987.5,1337.5,450,1655,389.294,389.722],'cm^-1')),
        HinderedRotor(inertia=(0.00111326,'amu*angstrom^2'), symmetry=1, barrier=(0.119627,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.00111174,'amu*angstrom^2'), symmetry=1, barrier=(0.119627,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.169593,'amu*angstrom^2'), symmetry=1, barrier=(18.1826,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.16914,'amu*angstrom^2'), symmetry=1, barrier=(18.1879,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.986737,'amu*angstrom^2'), symmetry=1, barrier=(105.511,'kJ/mol'), semiclassical=False),
    ],
    spinMultiplicity = 4,
    opticalIsomers = 1,
    molecularWeight = (97.1351,'amu'),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[0.332022,0.069614,-4.07632e-05,-5.28603e-09,9.77441e-12,31778.8,30.906], Tmin=(100,'K'), Tmax=(940.896,'K')), NASAPolynomial(coeffs=[18.208,0.0208187,-6.33558e-06,1.04518e-09,-7.14911e-14,27210.9,-60.6442], Tmin=(940.896,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(263.042,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(361.68,'J/(mol*K)'), comment="""Thermo group additivity estimation: group(O2s-(Cds-Cd)H) + group(Cs-CsCsHH) + group(Cs-(Cds-Cds)CsHH) + group(Cs-CsHHH) + group(Cs-(Cds-Cds)HHH) + group(Cds-CdsCsCs) + group(Cds-CdsOsH) + radical(RCCJC) + radical(RCCJ) + radical(Allyl_P)"""),
)

species(
    label = '[CH2]C([CH2])C(O)[C]=C(14106)',
    structure = SMILES('[CH2]C([CH2])C(O)[C]=C'),
    E0 = (402.076,'kJ/mol'),
    modes = [
        HarmonicOscillator(frequencies=([3000,3033.33,3066.67,3100,415,465,780,850,1435,1475,900,1100,2950,3100,1380,975,1025,1650,3615,1277.5,1000,1685,370,1380,1383.33,1386.67,1390,370,373.333,376.667,380,2800,3000,430,440,435.848,2021.89],'cm^-1')),
        HinderedRotor(inertia=(0.06962,'amu*angstrom^2'), symmetry=1, barrier=(9.25048,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.0681715,'amu*angstrom^2'), symmetry=1, barrier=(9.19838,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.0672189,'amu*angstrom^2'), symmetry=1, barrier=(9.26496,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.000894746,'amu*angstrom^2'), symmetry=1, barrier=(0.119627,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.522598,'amu*angstrom^2'), symmetry=1, barrier=(72.1474,'kJ/mol'), semiclassical=False),
    ],
    spinMultiplicity = 4,
    opticalIsomers = 1,
    molecularWeight = (97.1351,'amu'),
    collisionModel = TransportData(shapeIndex=2, epsilon=(3818.2,'J/mol'), sigma=(6.62498,'angstroms'), dipoleMoment=(0,'C*m'), polarizability=(0,'angstroms^3'), rotrelaxcollnum=0, comment="""Epsilon & sigma estimated with Tc=596.39 K, Pc=29.8 bar (from Joback method)"""),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[0.471104,0.0723196,-6.81821e-05,3.55797e-08,-7.47937e-12,48490.1,33.5628], Tmin=(100,'K'), Tmax=(1153.67,'K')), NASAPolynomial(coeffs=[13.5471,0.0269817,-9.2323e-06,1.51378e-09,-9.71143e-14,45473.1,-31.39], Tmin=(1153.67,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(402.076,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(361.68,'J/(mol*K)'), comment="""Thermo group additivity estimation: group(O2s-CsH) + group(Cs-CsCsCsH) + group(Cs-(Cds-Cds)CsOsH) + group(Cs-CsHHH) + group(Cs-CsHHH) + group(Cds-CdsCsH) + group(Cds-CdsHH) + radical(Isobutyl) + radical(Isobutyl) + radical(Cds_S)"""),
)

species(
    label = '[CH2]C1CC(O)C1=C(28530)',
    structure = SMILES('[CH2]C1CC(O)C1=C'),
    E0 = (100.877,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
    molecularWeight = (97.1351,'amu'),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[10.4365,-0.0112529,0.000128602,-1.24682e-07,2.88304e-11,11786.9,-15.0057], Tmin=(100,'K'), Tmax=(1711.31,'K')), NASAPolynomial(coeffs=[73.853,0.0284899,-7.09949e-05,1.72591e-08,-1.28186e-12,-37442.8,-435.445], Tmin=(1711.31,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(100.877,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(374.151,'J/(mol*K)'), comment="""Thermo group additivity estimation: group(O2s-CsH) + group(Cs-(Cds-Cds)CsCsH) + group(Cs-(Cds-Cds)CsOsH) + group(Cs-CsCsHH) + group(Cs-CsHHH) + group(Cds-CdsCsCs) + group(Cds-CdsHH) + ring(methylenecyclobutane) + radical(Isobutyl)"""),
)

species(
    label = 'C=[C]C(O)CC=C(27749)',
    structure = SMILES('C=[C]C(O)CC=C'),
    E0 = (128.735,'kJ/mol'),
    modes = [
        HarmonicOscillator(frequencies=([3615,1277.5,1000,2950,3000,3050,3100,1330,1430,900,1050,1000,1050,1600,1700,2750,2850,1437.5,1250,1305,750,350,1380,1390,370,380,2900,435,1685,370,3010,987.5,1337.5,450,1655,349.768,350.354,352.146],'cm^-1')),
        HinderedRotor(inertia=(0.0013623,'amu*angstrom^2'), symmetry=1, barrier=(0.119627,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.14261,'amu*angstrom^2'), symmetry=1, barrier=(12.5301,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.143242,'amu*angstrom^2'), symmetry=1, barrier=(12.5323,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.144417,'amu*angstrom^2'), symmetry=1, barrier=(12.5557,'kJ/mol'), semiclassical=False),
    ],
    spinMultiplicity = 2,
    opticalIsomers = 1,
    molecularWeight = (97.1351,'amu'),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[0.655462,0.0644142,-4.55079e-05,1.47399e-08,-1.46028e-12,15611.3,30.7952], Tmin=(100,'K'), Tmax=(1209.51,'K')), NASAPolynomial(coeffs=[14.7712,0.0270543,-1.07367e-05,1.94702e-09,-1.33206e-13,11514.7,-42.809], Tmin=(1209.51,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(128.735,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(365.837,'J/(mol*K)'), comment="""Thermo group additivity estimation: group(O2s-CsH) + group(Cs-(Cds-Cds)CsOsH) + group(Cs-(Cds-Cds)CsHH) + group(Cds-CdsCsH) + group(Cds-CdsCsH) + group(Cds-CdsHH) + group(Cds-CdsHH) + radical(Cds_S)"""),
)

species(
    label = '[C]=C(584)',
    structure = SMILES('[C]=C'),
    E0 = (600.251,'kJ/mol'),
    modes = [
        HarmonicOscillator(frequencies=([2950,3100,1380,975,1025,1650],'cm^-1')),
    ],
    spinMultiplicity = 3,
    opticalIsomers = 1,
    molecularWeight = (26.0373,'amu'),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[3.94093,-0.00117598,1.80376e-05,-2.01208e-08,6.96659e-12,72197.9,5.25681], Tmin=(100,'K'), Tmax=(976.125,'K')), NASAPolynomial(coeffs=[3.93016,0.00536132,-1.98619e-06,3.69549e-10,-2.66221e-14,71890.7,3.724], Tmin=(976.125,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(600.251,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(83.1447,'J/(mol*K)'), comment="""Thermo group additivity estimation: group(Cds-CdsHH) + group(Cds-CdsHH) + radical(CdCdJ2_triplet)"""),
)

species(
    label = '[CH2][CH]C[CH]O(747)',
    structure = SMILES('[CH2][CH]C[CH]O'),
    E0 = (275.859,'kJ/mol'),
    modes = [
        HarmonicOscillator(frequencies=([3000,3100,440,815,1455,1000,2750,2850,1437.5,1250,1305,750,350,3615,1277.5,1000,3000,3050,390,425,1340,1360,335,370,1553.81,1556.01],'cm^-1')),
        HinderedRotor(inertia=(0.181251,'amu*angstrom^2'), symmetry=1, barrier=(4.16731,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.181288,'amu*angstrom^2'), symmetry=1, barrier=(4.16816,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.180401,'amu*angstrom^2'), symmetry=1, barrier=(4.14776,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.181288,'amu*angstrom^2'), symmetry=1, barrier=(4.16817,'kJ/mol'), semiclassical=False),
    ],
    spinMultiplicity = 4,
    opticalIsomers = 1,
    molecularWeight = (71.0978,'amu'),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[1.72539,0.0558307,-7.95632e-05,7.03156e-08,-2.46489e-11,33254.5,24.5429], Tmin=(100,'K'), Tmax=(863.575,'K')), NASAPolynomial(coeffs=[4.59712,0.0300846,-1.32276e-05,2.41864e-09,-1.62389e-13,33222.5,13.7964], Tmin=(863.575,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(275.859,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(266.063,'J/(mol*K)'), comment="""Thermo group additivity estimation: group(O2s-CsH) + group(Cs-CsCsHH) + group(Cs-CsCsHH) + group(Cs-CsOsHH) + group(Cs-CsHHH) + radical(RCCJ) + radical(CCsJOH) + radical(RCCJC)"""),
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
    label = '[CH2]C(O)[C]=C(5788)',
    structure = SMILES('[CH2]C(O)[C]=C'),
    E0 = (260.102,'kJ/mol'),
    modes = [
        HarmonicOscillator(frequencies=([3000,3100,440,815,1455,1000,2950,3100,1380,975,1025,1650,3615,1277.5,1000,1685,370,1380,1390,370,380,2900,435,321.366],'cm^-1')),
        HinderedRotor(inertia=(0.0016324,'amu*angstrom^2'), symmetry=1, barrier=(0.119627,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.166096,'amu*angstrom^2'), symmetry=1, barrier=(12.1709,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.166087,'amu*angstrom^2'), symmetry=1, barrier=(12.171,'kJ/mol'), semiclassical=False),
    ],
    spinMultiplicity = 3,
    opticalIsomers = 1,
    molecularWeight = (70.0898,'amu'),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[1.6858,0.0455507,-3.94521e-05,1.7704e-08,-3.16924e-12,31370.9,22.7623], Tmin=(100,'K'), Tmax=(1345.51,'K')), NASAPolynomial(coeffs=[11.8503,0.0153332,-5.76502e-06,1.01282e-09,-6.79554e-14,28635.6,-29.2919], Tmin=(1345.51,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(260.102,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(245.277,'J/(mol*K)'), comment="""Thermo group additivity estimation: group(O2s-CsH) + group(Cs-(Cds-Cds)CsOsH) + group(Cs-CsHHH) + group(Cds-CdsCsH) + group(Cds-CdsHH) + radical(Cds_S) + radical(CJCO)"""),
)

species(
    label = '[CH2][C]CC(O)[C]=C(28551)',
    structure = SMILES('[CH2][C]CC(O)[C]=C'),
    E0 = (654.414,'kJ/mol'),
    modes = [
        HarmonicOscillator(frequencies=([3615,1277.5,1000,2950,3100,1380,975,1025,1650,2750,2850,1437.5,1250,1305,750,350,3000,3100,440,815,1455,1000,1685,370,1380,1390,370,380,2900,435,200,800,1200,1600],'cm^-1')),
        HinderedRotor(inertia=(0.156089,'amu*angstrom^2'), symmetry=1, barrier=(3.5888,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.156089,'amu*angstrom^2'), symmetry=1, barrier=(3.5888,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.156089,'amu*angstrom^2'), symmetry=1, barrier=(3.5888,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.156089,'amu*angstrom^2'), symmetry=1, barrier=(3.5888,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.156089,'amu*angstrom^2'), symmetry=1, barrier=(3.5888,'kJ/mol'), semiclassical=False),
    ],
    spinMultiplicity = 5,
    opticalIsomers = 1,
    molecularWeight = (96.1271,'amu'),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[0.768519,0.0709524,-6.79672e-05,3.48734e-08,-7.3163e-12,78824.3,31.6016], Tmin=(100,'K'), Tmax=(1136.64,'K')), NASAPolynomial(coeffs=[12.7059,0.0289429,-1.25279e-05,2.35688e-09,-1.64387e-13,76110.6,-27.5184], Tmin=(1136.64,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(654.414,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(336.736,'J/(mol*K)'), comment="""Thermo group additivity estimation: group(O2s-CsH) + group(Cs-CsCsHH) + group(Cs-(Cds-Cds)CsOsH) + group(Cs-CsCsHH) + group(Cs-CsHHH) + group(Cds-CdsCsH) + group(Cds-CdsHH) + radical(CCJ2_triplet) + radical(RCCJ) + radical(Cds_S)"""),
)

species(
    label = '[CH][CH]CC(O)[C]=C(28552)',
    structure = SMILES('[CH][CH]CC(O)[C]=C'),
    E0 = (643.614,'kJ/mol'),
    modes = [
        HarmonicOscillator(frequencies=([3615,1277.5,1000,2950,3100,1380,975,1025,1650,3025,407.5,1350,352.5,2750,2850,1437.5,1250,1305,750,350,1380,1390,370,380,2900,435,1685,370,200,800,1000,1200,1400,1600],'cm^-1')),
        HinderedRotor(inertia=(0.156089,'amu*angstrom^2'), symmetry=1, barrier=(3.5888,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.156089,'amu*angstrom^2'), symmetry=1, barrier=(3.5888,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.156089,'amu*angstrom^2'), symmetry=1, barrier=(3.5888,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.156089,'amu*angstrom^2'), symmetry=1, barrier=(3.5888,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.156089,'amu*angstrom^2'), symmetry=1, barrier=(3.5888,'kJ/mol'), semiclassical=False),
    ],
    spinMultiplicity = 5,
    opticalIsomers = 1,
    molecularWeight = (96.1271,'amu'),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[1.02136,0.0692167,-7.22043e-05,4.44747e-08,-1.1614e-11,77513,32.7516], Tmin=(100,'K'), Tmax=(910.378,'K')), NASAPolynomial(coeffs=[8.85864,0.0347813,-1.5466e-05,2.9253e-09,-2.03966e-13,76086,-4.32286], Tmin=(910.378,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(643.614,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(336.736,'J/(mol*K)'), comment="""Thermo group additivity estimation: group(O2s-CsH) + group(Cs-CsCsHH) + group(Cs-(Cds-Cds)CsOsH) + group(Cs-CsCsHH) + group(Cs-CsHHH) + group(Cds-CdsCsH) + group(Cds-CdsHH) + radical(RCCJC) + radical(Cds_S) + radical(CCJ2_triplet)"""),
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
    E0 = (400.645,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS2',
    E0 = (509.64,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS3',
    E0 = (418.457,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS4',
    E0 = (553.709,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS5',
    E0 = (480.377,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS6',
    E0 = (400.645,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS7',
    E0 = (528.364,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS8',
    E0 = (583.227,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS9',
    E0 = (543.105,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS10',
    E0 = (515.337,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS11',
    E0 = (542.622,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS12',
    E0 = (542.622,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS13',
    E0 = (516.781,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS14',
    E0 = (544.976,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS15',
    E0 = (563.954,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS16',
    E0 = (531.771,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS17',
    E0 = (486.336,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS18',
    E0 = (485.654,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS19',
    E0 = (574.772,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS20',
    E0 = (741.205,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS21',
    E0 = (849.906,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS22',
    E0 = (670.965,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS23',
    E0 = (709.952,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS24',
    E0 = (812.352,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS25',
    E0 = (859.546,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS26',
    E0 = (478.893,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS27',
    E0 = (489.614,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS28',
    E0 = (495.12,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS29',
    E0 = (559.394,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS30',
    E0 = (408.93,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS31',
    E0 = (400.645,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS32',
    E0 = (910.426,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS33',
    E0 = (851.16,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS34',
    E0 = (866.219,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS35',
    E0 = (855.419,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

reaction(
    label = 'reaction1',
    reactants = ['[CH2][CH]CC(O)[C]=C(14200)'],
    products = ['[CH2]C=C(87)', 'C=C=CO(12571)'],
    transitionState = 'TS1',
    kinetics = Arrhenius(A=(5e+12,'s^-1'), n=0, Ea=(0,'kJ/mol'), T0=(1,'K'), Tmin=(300,'K'), Tmax=(1500,'K'), comment="""Exact match found for rate rule [RJJ]
Euclidian distance = 0
family: 1,4_Linear_birad_scission"""),
)

reaction(
    label = 'reaction2',
    reactants = ['H(8)', '[CH2][CH]CC(O)=C=C(28542)'],
    products = ['[CH2][CH]CC(O)[C]=C(14200)'],
    transitionState = 'TS2',
    kinetics = Arrhenius(A=(169.619,'m^3/(mol*s)'), n=1.605, Ea=(12.4249,'kJ/mol'), T0=(1,'K'), comment="""Estimated using an average for rate rule [Cds_Ca;HJ]
Euclidian distance = 0
family: R_Addition_MultipleBond"""),
)

reaction(
    label = 'reaction3',
    reactants = ['H(8)', 'C=[C]C(O)[CH]C=C(28543)'],
    products = ['[CH2][CH]CC(O)[C]=C(14200)'],
    transitionState = 'TS3',
    kinetics = Arrhenius(A=(8.22e+08,'cm^3/(mol*s)'), n=1.533, Ea=(7.77387,'kJ/mol'), T0=(1,'K'), comment="""From training reaction 192 used for Cd_R;HJ
Exact match found for rate rule [Cd_R;HJ]
Euclidian distance = 0
family: R_Addition_MultipleBond"""),
)

reaction(
    label = 'reaction4',
    reactants = ['H(8)', 'C#CC(O)C[CH][CH2](22868)'],
    products = ['[CH2][CH]CC(O)[C]=C(14200)'],
    transitionState = 'TS4',
    kinetics = Arrhenius(A=(1.255e+11,'cm^3/(mol*s)'), n=1.005, Ea=(13.1503,'kJ/mol'), T0=(1,'K'), comment="""From training reaction 138 used for Ct-H_Ct-Cs;HJ
Exact match found for rate rule [Ct-H_Ct-Cs;HJ]
Euclidian distance = 0
family: R_Addition_MultipleBond"""),
)

reaction(
    label = 'reaction5',
    reactants = ['[CH2][CH][CH2](6136)', 'C=C=CO(12571)'],
    products = ['[CH2][CH]CC(O)[C]=C(14200)'],
    transitionState = 'TS5',
    kinetics = Arrhenius(A=(0.00803594,'m^3/(mol*s)'), n=2.41733, Ea=(22.1495,'kJ/mol'), T0=(1,'K'), comment="""Estimated using an average for rate rule [Cds_Ca;CJ]
Euclidian distance = 0
Multiplied by reaction path degeneracy 2.0
family: R_Addition_MultipleBond"""),
)

reaction(
    label = 'reaction6',
    reactants = ['[CH2]C=C(87)', '[CH2][C]=CO(18753)'],
    products = ['[CH2][CH]CC(O)[C]=C(14200)'],
    transitionState = 'TS6',
    kinetics = Arrhenius(A=(0.246938,'m^3/(mol*s)'), n=2.00579, Ea=(56.3505,'kJ/mol'), T0=(1,'K'), comment="""Estimated using an average for rate rule [Cd_R;CJ]
Euclidian distance = 0
Multiplied by reaction path degeneracy 2.0
family: R_Addition_MultipleBond
Ea raised from 51.6 to 56.4 kJ/mol to match endothermicity of reaction."""),
)

reaction(
    label = 'reaction7',
    reactants = ['OH(D)(132)', '[CH2][CH]CC=C=C(17733)'],
    products = ['[CH2][CH]CC(O)[C]=C(14200)'],
    transitionState = 'TS7',
    kinetics = Arrhenius(A=(986500,'cm^3/(mol*s)'), n=2.037, Ea=(0,'kJ/mol'), T0=(1,'K'), comment="""Estimated using template [Cds_Ca;OJ_pri] for rate rule [Cds-CsH_Ca;OJ_pri]
Euclidian distance = 1.0
family: R_Addition_MultipleBond
Ea raised from -6.0 to 0 kJ/mol."""),
)

reaction(
    label = 'reaction8',
    reactants = ['[CH2][CH]CC(O)[C]=C(14200)'],
    products = ['[CH2][CH]CC(O)=C[CH2](14197)'],
    transitionState = 'TS8',
    kinetics = Arrhenius(A=(3.677e+10,'s^-1'), n=0.839, Ea=(182.581,'kJ/mol'), T0=(1,'K'), comment="""Estimated using template [R2H_S;Cd_rad_out_Cd;Cs_H_out_noH] for rate rule [R2H_S;Cd_rad_out_Cd;Cs_H_out_NDMustO]
Euclidian distance = 2.0
family: intra_H_migration"""),
)

reaction(
    label = 'reaction9',
    reactants = ['[CH2]C[CH]C(O)[C]=C(28544)'],
    products = ['[CH2][CH]CC(O)[C]=C(14200)'],
    transitionState = 'TS9',
    kinetics = Arrhenius(A=(1.91383e+06,'s^-1'), n=2.02407, Ea=(137.004,'kJ/mol'), T0=(1,'K'), comment="""Estimated using an average for rate rule [R2H_S;C_rad_out_H/NonDeC;XH_out]
Euclidian distance = 0
Multiplied by reaction path degeneracy 2.0
family: intra_H_migration"""),
)

reaction(
    label = 'reaction10',
    reactants = ['[CH]=CC(O)C[CH][CH2](14203)'],
    products = ['[CH2][CH]CC(O)[C]=C(14200)'],
    transitionState = 'TS10',
    kinetics = Arrhenius(A=(1.08e+06,'s^-1'), n=1.99, Ea=(105.437,'kJ/mol'), T0=(1,'K'), comment="""From training reaction 17 used for R2H_D;Cd_rad_out_singleH;Cd_H_out_singleNd
Exact match found for rate rule [R2H_D;Cd_rad_out_singleH;Cd_H_out_singleNd]
Euclidian distance = 0
family: intra_H_migration"""),
)

reaction(
    label = 'reaction13',
    reactants = ['[CH2][CH]CC(O)[C]=C(14200)'],
    products = ['[CH2][CH]CC([O])C=C(12801)'],
    transitionState = 'TS11',
    kinetics = Arrhenius(A=(2.4115e+09,'s^-1'), n=1.00333, Ea=(141.977,'kJ/mol'), T0=(1,'K'), comment="""Estimated using template [R3H_SS_Cs;Cd_rad_out_Cd;XH_out] for rate rule [R3H_SS_Cs;Cd_rad_out_Cd;O_H_out]
Euclidian distance = 1.0
family: intra_H_migration"""),
)

reaction(
    label = 'reaction12',
    reactants = ['[CH2][CH]CC(O)[C]=C(14200)'],
    products = ['[CH2][CH][CH]C(O)C=C(14199)'],
    transitionState = 'TS12',
    kinetics = Arrhenius(A=(4.823e+09,'s^-1'), n=1.00333, Ea=(141.977,'kJ/mol'), T0=(1,'K'), comment="""Estimated using an average for rate rule [R3H_SS_Cs;Cd_rad_out_Cd;XH_out]
Euclidian distance = 0
Multiplied by reaction path degeneracy 2.0
family: intra_H_migration"""),
)

reaction(
    label = 'reaction13',
    reactants = ['[CH2][CH]CC(O)[C]=C(14200)'],
    products = ['[CH2][C]=C(O)CC[CH2](28545)'],
    transitionState = 'TS13',
    kinetics = Arrhenius(A=(285601,'s^-1'), n=2.01653, Ea=(116.136,'kJ/mol'), T0=(1,'K'), comment="""Estimated using an average for rate rule [R3H_SS_Cs;Y_rad_out;XH_out]
Euclidian distance = 0
family: intra_H_migration"""),
)

reaction(
    label = 'reaction14',
    reactants = ['[CH2][CH]CC(O)[C]=C(14200)'],
    products = ['C=[C]C(O)[CH][CH]C(27904)'],
    transitionState = 'TS14',
    kinetics = Arrhenius(A=(2.50921e+06,'s^-1'), n=1.8375, Ea=(144.331,'kJ/mol'), T0=(1,'K'), comment="""Estimated using template [R3Hall;C_rad_out_2H;Cs_H_out_H/NonDeC] for rate rule [R3HJ;C_rad_out_2H;Cs_H_out_H/(NonDeC/O)]
Euclidian distance = 1.41421356237
Multiplied by reaction path degeneracy 2.0
family: intra_H_migration"""),
)

reaction(
    label = 'reaction15',
    reactants = ['[CH2]CCC([O])[C]=C(14202)'],
    products = ['[CH2][CH]CC(O)[C]=C(14200)'],
    transitionState = 'TS15',
    kinetics = Arrhenius(A=(3.50344e+06,'s^-1'), n=1.80068, Ea=(127.394,'kJ/mol'), T0=(1,'K'), comment="""Estimated using an average for rate rule [R4H_SSS;O_rad_out;XH_out]
Euclidian distance = 0
Multiplied by reaction path degeneracy 2.0
family: intra_H_migration"""),
)

reaction(
    label = 'reaction16',
    reactants = ['[CH2][CH]CC(O)[C]=C(14200)'],
    products = ['[CH2][C]=C(O)C[CH]C(27905)'],
    transitionState = 'TS16',
    kinetics = Arrhenius(A=(434.916,'s^-1'), n=2.67444, Ea=(131.125,'kJ/mol'), T0=(1,'K'), comment="""Estimated using template [R4Hall;C_rad_out_2H;XH_out] for rate rule [R4HJ_1;C_rad_out_2H;XH_out]
Euclidian distance = 1.0
family: intra_H_migration"""),
)

reaction(
    label = 'reaction17',
    reactants = ['[CH]=[C]C(O)CC[CH2](28546)'],
    products = ['[CH2][CH]CC(O)[C]=C(14200)'],
    transitionState = 'TS17',
    kinetics = Arrhenius(A=(272000,'s^-1'), n=1.9199, Ea=(33.0402,'kJ/mol'), T0=(1,'K'), comment="""Estimated using template [R5Hall;Cd_rad_out_singleH;XH_out] for rate rule [R5HJ_1;Cd_rad_out_singleH;XH_out]
Euclidian distance = 1.0
Multiplied by reaction path degeneracy 2.0
family: intra_H_migration"""),
)

reaction(
    label = 'reaction18',
    reactants = ['C=[C]C([O])C[CH]C(13187)'],
    products = ['[CH2][CH]CC(O)[C]=C(14200)'],
    transitionState = 'TS18',
    kinetics = Arrhenius(A=(7.8e+08,'s^-1'), n=0.775, Ea=(59.894,'kJ/mol'), T0=(1,'K'), comment="""Estimated using template [R5Hall;O_rad_out;Cs_H_out_2H] for rate rule [R5HJ_3;O_rad_out;Cs_H_out_2H]
Euclidian distance = 1.0
Multiplied by reaction path degeneracy 3.0
family: intra_H_migration"""),
)

reaction(
    label = 'reaction19',
    reactants = ['[CH]=[C]C(O)C[CH]C(27906)'],
    products = ['[CH2][CH]CC(O)[C]=C(14200)'],
    transitionState = 'TS19',
    kinetics = Arrhenius(A=(22.7193,'s^-1'), n=3.21897, Ea=(132.277,'kJ/mol'), T0=(1,'K'), comment="""Estimated using template [RnH;Cd_rad_out_singleH;Cs_H_out_2H] for rate rule [R6Hall;Cd_rad_out_singleH;Cs_H_out_2H]
Euclidian distance = 1.0
Multiplied by reaction path degeneracy 3.0
family: intra_H_migration"""),
)

reaction(
    label = 'reaction20',
    reactants = ['OH(D)(132)', '[CH2][C]=CC[CH][CH2](17739)'],
    products = ['[CH2][CH]CC(O)[C]=C(14200)'],
    transitionState = 'TS20',
    kinetics = Arrhenius(A=(3.05166e+07,'m^3/(mol*s)'), n=0.045, Ea=(0.1046,'kJ/mol'), T0=(1,'K'), comment="""Estimated using an average for rate rule [O_pri_rad;Y_rad]
Euclidian distance = 0
family: R_Recombination"""),
)

reaction(
    label = 'reaction21',
    reactants = ['H(8)', '[CH2][CH]CC([O])[C]=C(14207)'],
    products = ['[CH2][CH]CC(O)[C]=C(14200)'],
    transitionState = 'TS21',
    kinetics = Arrhenius(A=(5.00518e+06,'m^3/(mol*s)'), n=0.282325, Ea=(7.09479,'kJ/mol'), T0=(1,'K'), comment="""Estimated using average of templates [Y_rad;O_rad/NonDe] + [H_rad;O_sec_rad] for rate rule [H_rad;O_rad/NonDe]
Euclidian distance = 1.0
family: R_Recombination"""),
)

reaction(
    label = 'reaction22',
    reactants = ['[CH2][CH][CH2](6136)', '[CH2][C]=CO(18753)'],
    products = ['[CH2][CH]CC(O)[C]=C(14200)'],
    transitionState = 'TS22',
    kinetics = Arrhenius(A=(3.9578e+07,'m^3/(mol*s)'), n=-0.126319, Ea=(0,'kJ/mol'), T0=(1,'K'), comment="""Estimated using an average for rate rule [Y_rad;Y_rad]
Euclidian distance = 0
Multiplied by reaction path degeneracy 2.0
family: R_Recombination
Ea raised from -15.6 to -15.6 kJ/mol.
Ea raised from -15.6 to 0 kJ/mol."""),
)

reaction(
    label = 'reaction23',
    reactants = ['H(8)', '[CH2][C]=C(O)C[CH][CH2](28547)'],
    products = ['[CH2][CH]CC(O)[C]=C(14200)'],
    transitionState = 'TS23',
    kinetics = Arrhenius(A=(4.34078e+06,'m^3/(mol*s)'), n=0.278577, Ea=(0,'kJ/mol'), T0=(1,'K'), comment="""Estimated using an average for rate rule [Y_rad;H_rad]
Euclidian distance = 0
family: R_Recombination
Ea raised from -1.4 to 0 kJ/mol."""),
)

reaction(
    label = 'reaction24',
    reactants = ['H(8)', '[CH2][CH][CH]C(O)[C]=C(28548)'],
    products = ['[CH2][CH]CC(O)[C]=C(14200)'],
    transitionState = 'TS24',
    kinetics = Arrhenius(A=(4.34078e+06,'m^3/(mol*s)'), n=0.278577, Ea=(0,'kJ/mol'), T0=(1,'K'), comment="""Estimated using an average for rate rule [Y_rad;H_rad]
Euclidian distance = 0
family: R_Recombination
Ea raised from -1.4 to 0 kJ/mol."""),
)

reaction(
    label = 'reaction25',
    reactants = ['H(8)', '[CH]=[C]C(O)C[CH][CH2](28549)'],
    products = ['[CH2][CH]CC(O)[C]=C(14200)'],
    transitionState = 'TS25',
    kinetics = Arrhenius(A=(4.34078e+06,'m^3/(mol*s)'), n=0.278577, Ea=(0,'kJ/mol'), T0=(1,'K'), comment="""Estimated using an average for rate rule [Y_rad;H_rad]
Euclidian distance = 0
family: R_Recombination
Ea raised from -1.4 to 0 kJ/mol."""),
)

reaction(
    label = 'reaction26',
    reactants = ['[CH2][CH]CC(O)[C]=C(14200)'],
    products = ['[CH2]CCC(O)=C=C(28550)'],
    transitionState = 'TS26',
    kinetics = Arrhenius(A=(2.00399e+09,'s^-1'), n=0.37, Ea=(78.2471,'kJ/mol'), T0=(1,'K'), comment="""Estimated using average of templates [R3;Y_rad;XH_Rrad_De] + [R3radExo;Y_rad;XH_Rrad] for rate rule [R3radExo;Y_rad;XH_Rrad_De]
Euclidian distance = 1.0
family: Intra_Disproportionation"""),
)

reaction(
    label = 'reaction27',
    reactants = ['[CH2][CH]CC(O)[C]=C(14200)'],
    products = ['C=C[CH]C(O)C=C(14211)'],
    transitionState = 'TS27',
    kinetics = Arrhenius(A=(5.2748e+09,'s^-1'), n=0.37, Ea=(88.9686,'kJ/mol'), T0=(1,'K'), comment="""Estimated using average of templates [R3;Y_rad_De;XH_Rrad] + [R3radExo;Y_rad;XH_Rrad] for rate rule [R3radExo;Y_rad_De;XH_Rrad]
Euclidian distance = 1.0
Multiplied by reaction path degeneracy 2.0
family: Intra_Disproportionation"""),
)

reaction(
    label = 'reaction28',
    reactants = ['[CH2][CH]CC(O)[C]=C(14200)'],
    products = ['[CH2][CH]CC([CH2])=CO(14139)'],
    transitionState = 'TS28',
    kinetics = Arrhenius(A=(8.66e+11,'s^-1'), n=0.438, Ea=(94.4747,'kJ/mol'), T0=(1,'K'), comment="""Estimated using an average for rate rule [cCs(-HR!H)CJ;CdsJ;C]
Euclidian distance = 0
family: 1,2_shiftC"""),
)

reaction(
    label = 'reaction26',
    reactants = ['[CH2]C([CH2])C(O)[C]=C(14106)'],
    products = ['[CH2][CH]CC(O)[C]=C(14200)'],
    transitionState = 'TS29',
    kinetics = Arrhenius(A=(2.66e+08,'s^-1'), n=1.36, Ea=(157.318,'kJ/mol'), T0=(1,'K'), comment="""Estimated using an average for rate rule [cCsCJ;CsJ-HH;C]
Euclidian distance = 0
Multiplied by reaction path degeneracy 2.0
family: 1,2_shiftC"""),
)

reaction(
    label = 'reaction30',
    reactants = ['[CH2][CH]CC(O)[C]=C(14200)'],
    products = ['[CH2]C1CC(O)C1=C(28530)'],
    transitionState = 'TS30',
    kinetics = Arrhenius(A=(1.62e+12,'s^-1'), n=-0.305, Ea=(8.28432,'kJ/mol'), T0=(1,'K'), comment="""Estimated using an average for rate rule [R4_SSS;Y_rad_out;Ypri_rad_out]
Euclidian distance = 0
family: Birad_recombination"""),
)

reaction(
    label = 'reaction31',
    reactants = ['[CH2][CH]CC(O)[C]=C(14200)'],
    products = ['C=[C]C(O)CC=C(27749)'],
    transitionState = 'TS31',
    kinetics = Arrhenius(A=(6.31e+07,'s^-1'), n=0, Ea=(0,'kJ/mol'), T0=(1,'K'), Tmin=(300,'K'), Tmax=(1500,'K'), comment="""From training reaction 2 used for Y_12_10
Exact match found for rate rule [Y_12_10]
Euclidian distance = 0
family: 1,2-Birad_to_alkene"""),
)

reaction(
    label = 'reaction32',
    reactants = ['[C]=C(584)', '[CH2][CH]C[CH]O(747)'],
    products = ['[CH2][CH]CC(O)[C]=C(14200)'],
    transitionState = 'TS32',
    kinetics = Arrhenius(A=(1.14854e+06,'m^3/(mol*s)'), n=0.575199, Ea=(34.3157,'kJ/mol'), T0=(1,'K'), comment="""Estimated using template [Y_rad;Birad] for rate rule [C_rad/H/CsO;Birad]
Euclidian distance = 4.0
family: Birad_R_Recombination"""),
)

reaction(
    label = 'reaction33',
    reactants = ['[CH][CH2](721)', '[CH2]C(O)[C]=C(5788)'],
    products = ['[CH2][CH]CC(O)[C]=C(14200)'],
    transitionState = 'TS33',
    kinetics = Arrhenius(A=(1.14854e+06,'m^3/(mol*s)'), n=0.575199, Ea=(34.3157,'kJ/mol'), T0=(1,'K'), comment="""Estimated using template [Y_rad;Birad] for rate rule [C_rad/H2/Cs;Birad]
Euclidian distance = 3.0
family: Birad_R_Recombination"""),
)

reaction(
    label = 'reaction34',
    reactants = ['H(8)', '[CH2][C]CC(O)[C]=C(28551)'],
    products = ['[CH2][CH]CC(O)[C]=C(14200)'],
    transitionState = 'TS34',
    kinetics = Arrhenius(A=(1e+07,'m^3/(mol*s)'), n=0, Ea=(0,'kJ/mol'), T0=(1,'K'), comment="""Estimated using an average for rate rule [H_rad;Birad]
Euclidian distance = 0
family: Birad_R_Recombination"""),
)

reaction(
    label = 'reaction35',
    reactants = ['H(8)', '[CH][CH]CC(O)[C]=C(28552)'],
    products = ['[CH2][CH]CC(O)[C]=C(14200)'],
    transitionState = 'TS35',
    kinetics = Arrhenius(A=(1e+07,'m^3/(mol*s)'), n=0, Ea=(0,'kJ/mol'), T0=(1,'K'), comment="""Estimated using an average for rate rule [H_rad;Birad]
Euclidian distance = 0
family: Birad_R_Recombination"""),
)

network(
    label = '5137',
    isomers = [
        '[CH2][CH]CC(O)[C]=C(14200)',
    ],
    reactants = [
        ('[CH2]C=C(87)', 'C=C=CO(12571)'),
    ],
    bathGas = {
        'N2': 0.25,
        'Ne': 0.25,
        'He': 0.25,
        'Ar': 0.25,
    },
)

pressureDependence(
    label = '5137',
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

