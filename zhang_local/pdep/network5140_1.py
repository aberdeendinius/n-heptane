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
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[0.471104,0.0723196,-6.81821e-05,3.55797e-08,-7.47937e-12,48490.1,33.5628], Tmin=(100,'K'), Tmax=(1153.67,'K')), NASAPolynomial(coeffs=[13.5471,0.0269817,-9.2323e-06,1.51378e-09,-9.71143e-14,45473.1,-31.39], Tmin=(1153.67,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(402.076,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(361.68,'J/(mol*K)'), comment="""Thermo group additivity estimation: group(O2s-CsH) + group(Cs-CsCsCsH) + group(Cs-(Cds-Cds)CsOsH) + group(Cs-CsHHH) + group(Cs-CsHHH) + group(Cds-CdsCsH) + group(Cds-CdsHH) + radical(Isobutyl) + radical(Isobutyl) + radical(Cds_S)"""),
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
    label = '[CH2]C(=C)C(O)[C]=C(18301)',
    structure = SMILES('[CH2]C(=C)C(O)[C]=C'),
    E0 = (261.203,'kJ/mol'),
    modes = [
        HarmonicOscillator(frequencies=([3615,1277.5,1000,350,440,435,1725,2950,3000,3050,3100,1330,1430,900,1050,1000,1050,1600,1700,3000,3100,440,815,1455,1000,1685,370,1380,1390,370,380,2900,435,180,814.927],'cm^-1')),
        HinderedRotor(inertia=(0.00667524,'amu*angstrom^2'), symmetry=1, barrier=(3.12414,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.791632,'amu*angstrom^2'), symmetry=1, barrier=(18.2012,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.792961,'amu*angstrom^2'), symmetry=1, barrier=(18.2317,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.792018,'amu*angstrom^2'), symmetry=1, barrier=(18.21,'kJ/mol'), semiclassical=False),
    ],
    spinMultiplicity = 3,
    opticalIsomers = 1,
    molecularWeight = (96.1271,'amu'),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[0.679628,0.066847,-5.5719e-05,2.3892e-08,-4.12578e-12,31540,28.6411], Tmin=(100,'K'), Tmax=(1377.62,'K')), NASAPolynomial(coeffs=[15.0824,0.0250272,-1.01835e-05,1.85588e-09,-1.26783e-13,27571.7,-45.4579], Tmin=(1377.62,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(261.203,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(340.893,'J/(mol*K)'), comment="""Thermo group additivity estimation: group(O2s-CsH) + group(Cs-(Cds-Cds)(Cds-Cds)OsH) + group(Cs-(Cds-Cds)HHH) + group(Cds-CdsCsCs) + group(Cds-CdsCsH) + group(Cds-CdsHH) + group(Cds-CdsHH) + radical(Cds_S) + radical(Allyl_P)"""),
)

species(
    label = '[CH2]C([CH2])C(O)=C=C(28507)',
    structure = SMILES('[CH2]C([CH2])C(O)=C=C'),
    E0 = (287.911,'kJ/mol'),
    modes = [
        HarmonicOscillator(frequencies=([3615,1277.5,1000,540,610,2055,350,440,435,1725,2950,3100,1380,975,1025,1650,3000,3033.33,3066.67,3100,415,465,780,850,1435,1475,900,1100,1380,1390,370,380,2900,435,552.862],'cm^-1')),
        HinderedRotor(inertia=(0.41853,'amu*angstrom^2'), symmetry=1, barrier=(13.0122,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.0609898,'amu*angstrom^2'), symmetry=1, barrier=(13.0193,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.0607963,'amu*angstrom^2'), symmetry=1, barrier=(13.0102,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(3.43287,'amu*angstrom^2'), symmetry=1, barrier=(78.9284,'kJ/mol'), semiclassical=False),
    ],
    spinMultiplicity = 3,
    opticalIsomers = 1,
    molecularWeight = (96.1271,'amu'),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[-0.275902,0.0780671,-7.82861e-05,3.98465e-08,-7.67622e-12,34795.3,31.2116], Tmin=(100,'K'), Tmax=(1454.21,'K')), NASAPolynomial(coeffs=[18.7677,0.0151833,-2.58972e-06,1.78291e-10,-2.93763e-15,30367.1,-63.9758], Tmin=(1454.21,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(287.911,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(340.893,'J/(mol*K)'), comment="""Thermo group additivity estimation: group(O2s-(Cds-Cd)H) + group(Cs-(Cds-Cdd-Cd)CsCsH) + group(Cs-CsHHH) + group(Cs-CsHHH) + group(Cds-CdsCsOs) + group(Cds-CdsHH) + group(Cdd-CdsCds) + radical(Isobutyl) + radical(Isobutyl)"""),
)

species(
    label = 'C#CC(O)C([CH2])[CH2](22792)',
    structure = SMILES('C#CC(O)C([CH2])[CH2]'),
    E0 = (330.184,'kJ/mol'),
    modes = [
        HarmonicOscillator(frequencies=([3000,3033.33,3066.67,3100,415,465,780,850,1435,1475,900,1100,750,770,3400,2100,3615,1277.5,1000,2175,525,1380,1383.33,1386.67,1390,370,373.333,376.667,380,2800,3000,430,440,2895.72],'cm^-1')),
        HinderedRotor(inertia=(0.0700618,'amu*angstrom^2'), symmetry=1, barrier=(9.91025,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.547414,'amu*angstrom^2'), symmetry=1, barrier=(77.4382,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.0700641,'amu*angstrom^2'), symmetry=1, barrier=(9.91026,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.547391,'amu*angstrom^2'), symmetry=1, barrier=(77.4379,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.547504,'amu*angstrom^2'), symmetry=1, barrier=(77.4379,'kJ/mol'), semiclassical=False),
    ],
    spinMultiplicity = 3,
    opticalIsomers = 1,
    molecularWeight = (96.1271,'amu'),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[0.174555,0.074267,-7.55929e-05,4.09742e-08,-8.56814e-12,39858.2,31.7651], Tmin=(100,'K'), Tmax=(1308.95,'K')), NASAPolynomial(coeffs=[15.6607,0.0196125,-4.56071e-06,5.18015e-10,-2.41553e-14,36432,-44.7172], Tmin=(1308.95,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(330.184,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(336.736,'J/(mol*K)'), comment="""Thermo group additivity estimation: group(O2s-CsH) + group(Cs-CsCsCsH) + group(Cs-(Cds-Cds)CsOsH) + group(Cs-CsHHH) + group(Cs-CsHHH) + group(Ct-CtCs) + group(Ct-CtH) + radical(Isobutyl) + radical(Isobutyl)"""),
)

species(
    label = 'CH2(T)(28)',
    structure = SMILES('[CH2]'),
    E0 = (381.37,'kJ/mol'),
    modes = [
        HarmonicOscillator(frequencies=([1066.91,2790.99,3622.37],'cm^-1')),
    ],
    spinMultiplicity = 3,
    opticalIsomers = 1,
    molecularWeight = (14.0266,'amu'),
    collisionModel = TransportData(shapeIndex=2, epsilon=(1197.29,'J/mol'), sigma=(3.8,'angstroms'), dipoleMoment=(0,'C*m'), polarizability=(0,'angstroms^3'), rotrelaxcollnum=0.0, comment="""GRI-Mech"""),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[4.01192,-0.000154979,3.26298e-06,-2.40422e-09,5.69497e-13,45867.7,0.5332], Tmin=(100,'K'), Tmax=(1104.58,'K')), NASAPolynomial(coeffs=[3.14983,0.00296674,-9.76056e-07,1.54115e-10,-9.50338e-15,46058.1,4.77808], Tmin=(1104.58,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(381.37,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(58.2013,'J/(mol*K)'), label="""CH2(T)""", comment="""Thermo library: primaryThermoLibrary"""),
)

species(
    label = 'C=[C]C(O)C=C(27708)',
    structure = SMILES('C=[C]C(O)C=C'),
    E0 = (148.758,'kJ/mol'),
    modes = [
        HarmonicOscillator(frequencies=([1685,370,2950,3000,3050,3100,1330,1430,900,1050,1000,1050,1600,1700,3010,987.5,1337.5,450,1655,3615,1277.5,1000,1380,1390,370,380,2900,435,401.788,402.949],'cm^-1')),
        HinderedRotor(inertia=(0.102172,'amu*angstrom^2'), symmetry=1, barrier=(11.7321,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.102173,'amu*angstrom^2'), symmetry=1, barrier=(11.7215,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.101377,'amu*angstrom^2'), symmetry=1, barrier=(11.7041,'kJ/mol'), semiclassical=False),
    ],
    spinMultiplicity = 2,
    opticalIsomers = 1,
    molecularWeight = (83.1085,'amu'),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[1.43816,0.050265,-3.66185e-05,1.35366e-08,-2.02373e-12,17988.7,25.5802], Tmin=(100,'K'), Tmax=(1564.23,'K')), NASAPolynomial(coeffs=[12.9464,0.0208365,-8.39847e-06,1.50935e-09,-1.0151e-13,14388.3,-35.0892], Tmin=(1564.23,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(148.758,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(295.164,'J/(mol*K)'), comment="""Thermo group additivity estimation: group(O2s-CsH) + group(Cs-(Cds-Cds)(Cds-Cds)OsH) + group(Cds-CdsCsH) + group(Cds-CdsCsH) + group(Cds-CdsHH) + group(Cds-CdsHH) + radical(Cds_S)"""),
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
    label = '[CH2]C([CH2])C=C=C(17716)',
    structure = SMILES('[CH2]C([CH2])C=C=C'),
    E0 = (502.469,'kJ/mol'),
    modes = [
        HarmonicOscillator(frequencies=([3000,3033.33,3066.67,3100,415,465,780,850,1435,1475,900,1100,2950,3100,1380,975,1025,1650,540,610,2055,1380,1390,370,380,2900,435,3010,987.5,1337.5,450,1655,241.346],'cm^-1')),
        HinderedRotor(inertia=(1.80686,'amu*angstrom^2'), symmetry=1, barrier=(77.6744,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.288944,'amu*angstrom^2'), symmetry=1, barrier=(12.7063,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.0121146,'amu*angstrom^2'), symmetry=1, barrier=(77.6103,'kJ/mol'), semiclassical=False),
    ],
    spinMultiplicity = 3,
    opticalIsomers = 1,
    molecularWeight = (80.1277,'amu'),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[1.41477,0.0511993,-2.71882e-05,-2.1312e-09,5.51936e-12,60531.4,25.0871], Tmin=(100,'K'), Tmax=(922.534,'K')), NASAPolynomial(coeffs=[10.8754,0.0248363,-8.15522e-06,1.33679e-09,-8.75163e-14,58162.1,-23.1725], Tmin=(922.534,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(502.469,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(320.107,'J/(mol*K)'), comment="""Thermo group additivity estimation: group(Cs-(Cds-Cdd-Cd)CsCsH) + group(Cs-CsHHH) + group(Cs-CsHHH) + group(Cds-CdsCsH) + group(Cds-CdsHH) + group(Cdd-CdsCds) + radical(Isobutyl) + radical(Isobutyl)"""),
)

species(
    label = '[CH2][C](C)C(O)[C]=C(27854)',
    structure = SMILES('[CH2][C](C)C(O)[C]=C'),
    E0 = (349.568,'kJ/mol'),
    modes = [
        HarmonicOscillator(frequencies=([360,370,350,3615,1277.5,1000,2750,2800,2850,1350,1500,750,1050,1375,1000,2950,3100,1380,975,1025,1650,3000,3100,440,815,1455,1000,1685,370,1380,1390,370,380,2900,435,180,855.861],'cm^-1')),
        HinderedRotor(inertia=(0.115058,'amu*angstrom^2'), symmetry=1, barrier=(2.6454,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.000560641,'amu*angstrom^2'), symmetry=1, barrier=(2.64506,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.0347065,'amu*angstrom^2'), symmetry=1, barrier=(18.0407,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.784651,'amu*angstrom^2'), symmetry=1, barrier=(18.0407,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(3.47936,'amu*angstrom^2'), symmetry=1, barrier=(79.9974,'kJ/mol'), semiclassical=False),
    ],
    spinMultiplicity = 4,
    opticalIsomers = 1,
    molecularWeight = (97.1351,'amu'),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[0.74394,0.0637272,-4.32883e-05,1.18145e-08,-2.73205e-14,42167.2,31.6925], Tmin=(100,'K'), Tmax=(1071.04,'K')), NASAPolynomial(coeffs=[13.3347,0.0285842,-1.07076e-05,1.89069e-09,-1.28209e-13,38788.8,-33.0954], Tmin=(1071.04,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(349.568,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(361.68,'J/(mol*K)'), comment="""Thermo group additivity estimation: group(O2s-CsH) + group(Cs-CsCsCsH) + group(Cs-(Cds-Cds)CsOsH) + group(Cs-CsHHH) + group(Cs-CsHHH) + group(Cds-CdsCsH) + group(Cds-CdsHH) + radical(CCJ(C)CO) + radical(Isobutyl) + radical(Cds_S)"""),
)

species(
    label = '[CH2]C=C(O)C([CH2])[CH2](14104)',
    structure = SMILES('[CH2]C=C(O)C([CH2])[CH2]'),
    E0 = (262.806,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
    molecularWeight = (97.1351,'amu'),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[-0.70468,0.0802958,-7.50109e-05,3.54396e-08,-6.29781e-12,31797.5,33.9931], Tmin=(100,'K'), Tmax=(1613.16,'K')), NASAPolynomial(coeffs=[19.5277,0.016236,-2.52776e-06,1.46493e-10,-9.61616e-16,27077.4,-67.6891], Tmin=(1613.16,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(262.806,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(361.68,'J/(mol*K)'), comment="""Thermo group additivity estimation: group(O2s-(Cds-Cd)H) + group(Cs-(Cds-Cds)CsCsH) + group(Cs-CsHHH) + group(Cs-CsHHH) + group(Cs-(Cds-Cds)HHH) + group(Cds-CdsCsOs) + group(Cds-CdsCsH) + radical(Allyl_P) + radical(Isobutyl) + radical(Isobutyl)"""),
)

species(
    label = '[CH]=CC(O)C([CH2])[CH2](14107)',
    structure = SMILES('[CH]=CC(O)C([CH2])[CH2]'),
    E0 = (411.33,'kJ/mol'),
    modes = [
        HarmonicOscillator(frequencies=([3000,3033.33,3066.67,3100,415,465,780,850,1435,1475,900,1100,3120,650,792.5,1650,3010,987.5,1337.5,450,1655,3615,1277.5,1000,1380,1383.33,1386.67,1390,370,373.333,376.667,380,2800,3000,430,440,334.815],'cm^-1')),
        HinderedRotor(inertia=(0.00150802,'amu*angstrom^2'), symmetry=1, barrier=(0.119627,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.00150612,'amu*angstrom^2'), symmetry=1, barrier=(0.119627,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.0718407,'amu*angstrom^2'), symmetry=1, barrier=(5.69378,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.0717336,'amu*angstrom^2'), symmetry=1, barrier=(5.70006,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(1.00079,'amu*angstrom^2'), symmetry=1, barrier=(79.5764,'kJ/mol'), semiclassical=False),
    ],
    spinMultiplicity = 4,
    opticalIsomers = 1,
    molecularWeight = (97.1351,'amu'),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[0.227234,0.0736022,-6.81359e-05,3.39386e-08,-6.70021e-12,49615.4,34.3695], Tmin=(100,'K'), Tmax=(1273.54,'K')), NASAPolynomial(coeffs=[15.6976,0.0234608,-7.25125e-06,1.11061e-09,-6.82007e-14,45800.8,-43.5129], Tmin=(1273.54,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(411.33,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(361.68,'J/(mol*K)'), comment="""Thermo group additivity estimation: group(O2s-CsH) + group(Cs-CsCsCsH) + group(Cs-(Cds-Cds)CsOsH) + group(Cs-CsHHH) + group(Cs-CsHHH) + group(Cds-CdsCsH) + group(Cds-CdsHH) + radical(Isobutyl) + radical(Isobutyl) + radical(Cds_P)"""),
)

species(
    label = '[CH2][C]=C(O)C([CH2])C(27855)',
    structure = SMILES('[CH2][C]=C(O)C([CH2])C'),
    E0 = (295.565,'kJ/mol'),
    modes = [
        HarmonicOscillator(frequencies=([3615,1277.5,1000,2750,2800,2850,1350,1500,750,1050,1375,1000,350,440,435,1725,3000,3033.33,3066.67,3100,415,465,780,850,1435,1475,900,1100,1685,370,1380,1390,370,380,2900,435,180],'cm^-1')),
        HinderedRotor(inertia=(0.0292719,'amu*angstrom^2'), symmetry=1, barrier=(13.5024,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.587269,'amu*angstrom^2'), symmetry=1, barrier=(13.5025,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.587268,'amu*angstrom^2'), symmetry=1, barrier=(13.5024,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.0292735,'amu*angstrom^2'), symmetry=1, barrier=(13.5025,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.244257,'amu*angstrom^2'), symmetry=1, barrier=(112.644,'kJ/mol'), semiclassical=False),
    ],
    spinMultiplicity = 4,
    opticalIsomers = 1,
    molecularWeight = (97.1351,'amu'),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[-0.416015,0.0792659,-7.37195e-05,3.48084e-08,-6.29029e-12,35722.5,32.3378], Tmin=(100,'K'), Tmax=(1518.87,'K')), NASAPolynomial(coeffs=[19.7487,0.0174946,-4.15632e-06,5.18663e-10,-2.79609e-14,30596.7,-70.0824], Tmin=(1518.87,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(295.565,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(361.68,'J/(mol*K)'), comment="""Thermo group additivity estimation: group(O2s-(Cds-Cd)H) + group(Cs-(Cds-Cds)CsCsH) + group(Cs-CsHHH) + group(Cs-CsHHH) + group(Cs-(Cds-Cds)HHH) + group(Cds-CdsCsOs) + group(Cds-CdsCsH) + radical(Isobutyl) + radical(Allyl_P) + radical(Cds_S)"""),
)

species(
    label = '[CH2]C([CH2])C([O])C=C(12805)',
    structure = SMILES('[CH2]C([CH2])C([O])C=C'),
    E0 = (394.595,'kJ/mol'),
    modes = [
        HarmonicOscillator(frequencies=([3000,3033.33,3066.67,3100,415,465,780,850,1435,1475,900,1100,2950,3100,1380,975,1025,1650,1380,1383.33,1386.67,1390,370,373.333,376.667,380,2800,3000,430,440,3010,987.5,1337.5,450,1655,180,691.491,691.503],'cm^-1')),
        HinderedRotor(inertia=(0.0371361,'amu*angstrom^2'), symmetry=1, barrier=(12.6002,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.218081,'amu*angstrom^2'), symmetry=1, barrier=(73.9979,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.188078,'amu*angstrom^2'), symmetry=1, barrier=(4.32429,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.0051556,'amu*angstrom^2'), symmetry=1, barrier=(12.6002,'kJ/mol'), semiclassical=False),
    ],
    spinMultiplicity = 4,
    opticalIsomers = 1,
    molecularWeight = (97.1351,'amu'),
    collisionModel = TransportData(shapeIndex=2, epsilon=(3818.2,'J/mol'), sigma=(6.62498,'angstroms'), dipoleMoment=(0,'C*m'), polarizability=(0,'angstroms^3'), rotrelaxcollnum=0, comment="""Epsilon & sigma estimated with Tc=596.39 K, Pc=29.8 bar (from Joback method)"""),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[0.762047,0.0625905,-3.34704e-05,-3.96873e-09,7.33263e-12,47583.2,31.8543], Tmin=(100,'K'), Tmax=(949.086,'K')), NASAPolynomial(coeffs=[14.1898,0.0265003,-8.83396e-06,1.48667e-09,-9.98439e-14,44111,-37.0896], Tmin=(949.086,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(394.595,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(365.837,'J/(mol*K)'), comment="""Thermo group additivity estimation: group(O2s-CsH) + group(Cs-CsCsCsH) + group(Cs-(Cds-Cds)CsOsH) + group(Cs-CsHHH) + group(Cs-CsHHH) + group(Cds-CdsCsH) + group(Cds-CdsHH) + radical(Isobutyl) + radical(Isobutyl) + radical(CC(C)OJ)"""),
)

species(
    label = '[CH2][C]([CH2])C(O)C=C(14105)',
    structure = SMILES('[CH2][C]([CH2])C(O)C=C'),
    E0 = (316.808,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
    molecularWeight = (97.1351,'amu'),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[0.877269,0.0598618,-2.77846e-05,-8.95353e-09,8.95226e-12,38223.6,31.8278], Tmin=(100,'K'), Tmax=(945.617,'K')), NASAPolynomial(coeffs=[13.7665,0.0266243,-8.82358e-06,1.48165e-09,-9.95166e-14,34834.4,-34.6661], Tmin=(945.617,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(316.808,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(361.68,'J/(mol*K)'), comment="""Thermo group additivity estimation: group(O2s-CsH) + group(Cs-CsCsCsH) + group(Cs-(Cds-Cds)CsOsH) + group(Cs-CsHHH) + group(Cs-CsHHH) + group(Cds-CdsCsH) + group(Cds-CdsHH) + radical(Isobutyl) + radical(CCJ(C)CO) + radical(Isobutyl)"""),
)

species(
    label = '[CH2]C(C)C([O])[C]=C(13089)',
    structure = SMILES('[CH2]C(C)C([O])[C]=C'),
    E0 = (427.354,'kJ/mol'),
    modes = [
        HarmonicOscillator(frequencies=([1685,370,2750,2800,2850,1350,1500,750,1050,1375,1000,2950,3100,1380,975,1025,1650,1380,1383.33,1386.67,1390,370,373.333,376.667,380,2800,3000,430,440,3000,3100,440,815,1455,1000,180,180,180],'cm^-1')),
        HinderedRotor(inertia=(0.000184687,'amu*angstrom^2'), symmetry=1, barrier=(2.09694,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.0569979,'amu*angstrom^2'), symmetry=1, barrier=(18.9089,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.822538,'amu*angstrom^2'), symmetry=1, barrier=(18.9118,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.30494,'amu*angstrom^2'), symmetry=1, barrier=(101.38,'kJ/mol'), semiclassical=False),
    ],
    spinMultiplicity = 4,
    opticalIsomers = 1,
    molecularWeight = (97.1351,'amu'),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[0.474143,0.0681588,-5.43555e-05,2.2942e-08,-3.9065e-12,51533.6,32.2814], Tmin=(100,'K'), Tmax=(1403.13,'K')), NASAPolynomial(coeffs=[15.2075,0.0261572,-9.45395e-06,1.60788e-09,-1.0532e-13,47399.1,-43.7889], Tmin=(1403.13,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(427.354,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(365.837,'J/(mol*K)'), comment="""Thermo group additivity estimation: group(O2s-CsH) + group(Cs-CsCsCsH) + group(Cs-(Cds-Cds)CsOsH) + group(Cs-CsHHH) + group(Cs-CsHHH) + group(Cds-CdsCsH) + group(Cds-CdsHH) + radical(CC(C)OJ) + radical(Isobutyl) + radical(Cds_S)"""),
)

species(
    label = '[CH]=[C]C(O)C([CH2])C(27856)',
    structure = SMILES('[CH]=[C]C(O)C([CH2])C'),
    E0 = (444.09,'kJ/mol'),
    modes = [
        HarmonicOscillator(frequencies=([3120,650,792.5,1650,3615,1277.5,1000,2750,2800,2850,1350,1500,750,1050,1375,1000,3000,3100,440,815,1455,1000,1685,370,1380,1383.33,1386.67,1390,370,373.333,376.667,380,2800,3000,430,440,347.826],'cm^-1')),
        HinderedRotor(inertia=(0.0213737,'amu*angstrom^2'), symmetry=1, barrier=(1.83775,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.110467,'amu*angstrom^2'), symmetry=1, barrier=(9.47094,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.0214501,'amu*angstrom^2'), symmetry=1, barrier=(1.83923,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.11037,'amu*angstrom^2'), symmetry=1, barrier=(9.47115,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.200392,'amu*angstrom^2'), symmetry=1, barrier=(17.1797,'kJ/mol'), semiclassical=False),
    ],
    spinMultiplicity = 4,
    opticalIsomers = 1,
    molecularWeight = (97.1351,'amu'),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[0.499355,0.072722,-6.72206e-05,3.36791e-08,-6.82516e-12,53541.3,32.7775], Tmin=(100,'K'), Tmax=(1187.84,'K')), NASAPolynomial(coeffs=[13.8425,0.0277893,-1.04795e-05,1.83339e-09,-1.22677e-13,50371.4,-33.8921], Tmin=(1187.84,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(444.09,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(361.68,'J/(mol*K)'), comment="""Thermo group additivity estimation: group(O2s-CsH) + group(Cs-CsCsCsH) + group(Cs-(Cds-Cds)CsOsH) + group(Cs-CsHHH) + group(Cs-CsHHH) + group(Cds-CdsCsH) + group(Cds-CdsHH) + radical(Cds_P) + radical(Cds_S) + radical(Isobutyl)"""),
)

species(
    label = '[CH2][C]=CC([CH2])[CH2](17719)',
    structure = SMILES('[CH2][C]=CC([CH2])[CH2]'),
    E0 = (715.206,'kJ/mol'),
    modes = [
        HarmonicOscillator(frequencies=([3000,3020,3040,3060,3080,3100,415,440,465,780,815,850,1435,1455,1475,900,1000,1100,1685,370,1380,1390,370,380,2900,435,3010,987.5,1337.5,450,1655,2680.74],'cm^-1')),
        HinderedRotor(inertia=(0.0520452,'amu*angstrom^2'), symmetry=1, barrier=(13.9389,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(3.40752,'amu*angstrom^2'), symmetry=1, barrier=(78.3456,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(3.40743,'amu*angstrom^2'), symmetry=1, barrier=(78.3436,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.0153603,'amu*angstrom^2'), symmetry=1, barrier=(78.3389,'kJ/mol'), semiclassical=False),
    ],
    spinMultiplicity = 5,
    opticalIsomers = 1,
    molecularWeight = (80.1277,'amu'),
    collisionModel = TransportData(shapeIndex=2, epsilon=(3311.72,'J/mol'), sigma=(5.95882,'angstroms'), dipoleMoment=(0,'C*m'), polarizability=(0,'angstroms^3'), rotrelaxcollnum=0, comment="""Epsilon & sigma estimated with Tc=517.28 K, Pc=35.52 bar (from Joback method)"""),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[1.02621,0.0573665,-4.72316e-05,2.18687e-08,-4.10204e-12,86133.3,28.1831], Tmin=(100,'K'), Tmax=(1338.17,'K')), NASAPolynomial(coeffs=[11.978,0.0231449,-7.20676e-06,1.09932e-09,-6.69163e-14,83335.2,-27.3465], Tmin=(1338.17,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(715.206,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(315.95,'J/(mol*K)'), comment="""Thermo group additivity estimation: group(Cs-(Cds-Cds)CsCsH) + group(Cs-CsHHH) + group(Cs-CsHHH) + group(Cs-(Cds-Cds)HHH) + group(Cds-CdsCsH) + group(Cds-CdsCsH) + radical(Allyl_P) + radical(Isobutyl) + radical(Isobutyl) + radical(Cds_S)"""),
)

species(
    label = '[CH2]C([CH2])C([O])[C]=C(14110)',
    structure = SMILES('[CH2]C([CH2])C([O])[C]=C'),
    E0 = (632.437,'kJ/mol'),
    modes = [
        HarmonicOscillator(frequencies=([1685,370,2950,3100,1380,975,1025,1650,1380,1383.33,1386.67,1390,370,373.333,376.667,380,2800,3000,430,440,3000,3033.33,3066.67,3100,415,465,780,850,1435,1475,900,1100,342.784,342.808,342.818],'cm^-1')),
        HinderedRotor(inertia=(0.972513,'amu*angstrom^2'), symmetry=1, barrier=(81.0947,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.00143421,'amu*angstrom^2'), symmetry=1, barrier=(0.119627,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(1.31736e-05,'amu*angstrom^2'), symmetry=1, barrier=(0.119627,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.97244,'amu*angstrom^2'), symmetry=1, barrier=(81.0956,'kJ/mol'), semiclassical=False),
    ],
    spinMultiplicity = 5,
    opticalIsomers = 1,
    molecularWeight = (96.1271,'amu'),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[0.5983,0.0687351,-6.36527e-05,3.22455e-08,-6.56092e-12,76192.2,32.915], Tmin=(100,'K'), Tmax=(1191.99,'K')), NASAPolynomial(coeffs=[13.6634,0.0248919,-8.48001e-06,1.38779e-09,-8.89867e-14,73077.5,-32.4109], Tmin=(1191.99,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(632.437,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(340.893,'J/(mol*K)'), comment="""Thermo group additivity estimation: group(O2s-CsH) + group(Cs-CsCsCsH) + group(Cs-(Cds-Cds)CsOsH) + group(Cs-CsHHH) + group(Cs-CsHHH) + group(Cds-CdsCsH) + group(Cds-CdsHH) + radical(Isobutyl) + radical(CC(C)OJ) + radical(Cds_S) + radical(Isobutyl)"""),
)

species(
    label = '[CH2][C]([CH2])C(O)[C]=C(28508)',
    structure = SMILES('[CH2][C]([CH2])C(O)[C]=C'),
    E0 = (554.65,'kJ/mol'),
    modes = [
        HarmonicOscillator(frequencies=([360,370,350,3615,1277.5,1000,2950,3100,1380,975,1025,1650,3000,3033.33,3066.67,3100,415,465,780,850,1435,1475,900,1100,1685,370,1380,1390,370,380,2900,435,346.217,2001.46],'cm^-1')),
        HinderedRotor(inertia=(0.025974,'amu*angstrom^2'), symmetry=1, barrier=(73.8344,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.00140639,'amu*angstrom^2'), symmetry=1, barrier=(0.119627,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.00140632,'amu*angstrom^2'), symmetry=1, barrier=(0.119627,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.157447,'amu*angstrom^2'), symmetry=1, barrier=(13.3923,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.868066,'amu*angstrom^2'), symmetry=1, barrier=(73.8344,'kJ/mol'), semiclassical=False),
    ],
    spinMultiplicity = 5,
    opticalIsomers = 1,
    molecularWeight = (96.1271,'amu'),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[0.671243,0.0664974,-5.96488e-05,2.93934e-08,-5.83186e-12,66834.5,33.0407], Tmin=(100,'K'), Tmax=(1220.26,'K')), NASAPolynomial(coeffs=[13.4116,0.0247345,-8.31166e-06,1.34618e-09,-8.56701e-14,63725.2,-30.9601], Tmin=(1220.26,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(554.65,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(336.736,'J/(mol*K)'), comment="""Thermo group additivity estimation: group(O2s-CsH) + group(Cs-CsCsCsH) + group(Cs-(Cds-Cds)CsOsH) + group(Cs-CsHHH) + group(Cs-CsHHH) + group(Cds-CdsCsH) + group(Cds-CdsHH) + radical(Isobutyl) + radical(CCJ(C)CO) + radical(Isobutyl) + radical(Cds_S)"""),
)

species(
    label = '[CH2][C]=C(O)C([CH2])[CH2](28509)',
    structure = SMILES('[CH2][C]=C(O)C([CH2])[CH2]'),
    E0 = (500.648,'kJ/mol'),
    modes = [
        HarmonicOscillator(frequencies=([3000,3020,3040,3060,3080,3100,415,440,465,780,815,850,1435,1455,1475,900,1000,1100,1685,370,3615,1277.5,1000,350,440,435,1725,1380,1390,370,380,2900,435,359.62],'cm^-1')),
        HinderedRotor(inertia=(0.00130354,'amu*angstrom^2'), symmetry=1, barrier=(0.119627,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.12371,'amu*angstrom^2'), symmetry=1, barrier=(11.3531,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.00130351,'amu*angstrom^2'), symmetry=1, barrier=(0.119627,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(1.10252,'amu*angstrom^2'), symmetry=1, barrier=(101.177,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(1.10249,'amu*angstrom^2'), symmetry=1, barrier=(101.177,'kJ/mol'), semiclassical=False),
    ],
    spinMultiplicity = 5,
    opticalIsomers = 1,
    molecularWeight = (96.1271,'amu'),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[-0.252523,0.0794453,-8.18898e-05,4.29407e-08,-8.54375e-12,60379.2,32.8251], Tmin=(100,'K'), Tmax=(1396.55,'K')), NASAPolynomial(coeffs=[18.6806,0.0155166,-2.80669e-06,2.15389e-10,-5.05052e-15,56036.9,-61.4535], Tmin=(1396.55,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(500.648,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(336.736,'J/(mol*K)'), comment="""Thermo group additivity estimation: group(O2s-(Cds-Cd)H) + group(Cs-(Cds-Cds)CsCsH) + group(Cs-CsHHH) + group(Cs-CsHHH) + group(Cs-(Cds-Cds)HHH) + group(Cds-CdsCsOs) + group(Cds-CdsCsH) + radical(Isobutyl) + radical(Cds_S) + radical(Allyl_P) + radical(Isobutyl)"""),
)

species(
    label = '[CH]=[C]C(O)C([CH2])[CH2](28510)',
    structure = SMILES('[CH]=[C]C(O)C([CH2])[CH2]'),
    E0 = (649.172,'kJ/mol'),
    modes = [
        HarmonicOscillator(frequencies=([3000,3033.33,3066.67,3100,415,465,780,850,1435,1475,900,1100,1685,370,3615,1277.5,1000,3120,650,792.5,1650,1380,1383.33,1386.67,1390,370,373.333,376.667,380,2800,3000,430,440,276.425],'cm^-1')),
        HinderedRotor(inertia=(0.206526,'amu*angstrom^2'), symmetry=1, barrier=(11.3603,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.00380707,'amu*angstrom^2'), symmetry=1, barrier=(11.3562,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.00221686,'amu*angstrom^2'), symmetry=1, barrier=(0.119627,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.0021901,'amu*angstrom^2'), symmetry=1, barrier=(0.119627,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(1.31247,'amu*angstrom^2'), symmetry=1, barrier=(71.1351,'kJ/mol'), semiclassical=False),
    ],
    spinMultiplicity = 5,
    opticalIsomers = 1,
    molecularWeight = (96.1271,'amu'),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[0.483296,0.0748253,-8.12875e-05,4.84001e-08,-1.14821e-11,78206.2,33.9228], Tmin=(100,'K'), Tmax=(1031.21,'K')), NASAPolynomial(coeffs=[13.3656,0.0248558,-8.60188e-06,1.40964e-09,-9.00521e-14,75549.3,-28.6228], Tmin=(1031.21,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(649.172,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(336.736,'J/(mol*K)'), comment="""Thermo group additivity estimation: group(O2s-CsH) + group(Cs-CsCsCsH) + group(Cs-(Cds-Cds)CsOsH) + group(Cs-CsHHH) + group(Cs-CsHHH) + group(Cds-CdsCsH) + group(Cds-CdsHH) + radical(Isobutyl) + radical(Cds_S) + radical(Cds_P) + radical(Isobutyl)"""),
)

species(
    label = 'C=[C]C(O)C(=C)C(27848)',
    structure = SMILES('C=[C]C(O)C(=C)C'),
    E0 = (109.703,'kJ/mol'),
    modes = [
        HarmonicOscillator(frequencies=([3615,1277.5,1000,2750,2800,2850,1350,1500,750,1050,1375,1000,350,440,435,1725,2950,3000,3050,3100,1330,1430,900,1050,1000,1050,1600,1700,1380,1390,370,380,2900,435,1685,370,283.309,284.086],'cm^-1')),
        HinderedRotor(inertia=(0.00209818,'amu*angstrom^2'), symmetry=1, barrier=(0.119627,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.204518,'amu*angstrom^2'), symmetry=1, barrier=(11.7446,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.208241,'amu*angstrom^2'), symmetry=1, barrier=(11.7478,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.203808,'amu*angstrom^2'), symmetry=1, barrier=(11.7333,'kJ/mol'), semiclassical=False),
    ],
    spinMultiplicity = 2,
    opticalIsomers = 1,
    molecularWeight = (97.1351,'amu'),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[0.991706,0.0657602,-5.21868e-05,2.21763e-08,-3.9458e-12,13302.7,27.9277], Tmin=(100,'K'), Tmax=(1296.94,'K')), NASAPolynomial(coeffs=[11.6329,0.0329407,-1.42287e-05,2.66456e-09,-1.84689e-13,10542.5,-26.1765], Tmin=(1296.94,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(109.703,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(365.837,'J/(mol*K)'), comment="""Thermo group additivity estimation: group(O2s-CsH) + group(Cs-(Cds-Cds)(Cds-Cds)OsH) + group(Cs-(Cds-Cds)HHH) + group(Cds-CdsCsCs) + group(Cds-CdsCsH) + group(Cds-CdsHH) + group(Cds-CdsHH) + radical(Cds_S)"""),
)

species(
    label = '[CH2]C(=C)C(O)C=C(14114)',
    structure = SMILES('[CH2]C(=C)C(O)C=C'),
    E0 = (23.3607,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
    molecularWeight = (97.1351,'amu'),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[0.657256,0.0629412,-3.36099e-05,-1.44424e-09,4.93388e-12,2939.09,28.2456], Tmin=(100,'K'), Tmax=(1052.98,'K')), NASAPolynomial(coeffs=[15.3657,0.0269594,-1.06895e-05,1.98496e-09,-1.39801e-13,-1261.23,-48.71], Tmin=(1052.98,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(23.3607,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(365.837,'J/(mol*K)'), comment="""Thermo group additivity estimation: group(O2s-CsH) + group(Cs-(Cds-Cds)(Cds-Cds)OsH) + group(Cs-(Cds-Cds)HHH) + group(Cds-CdsCsCs) + group(Cds-CdsCsH) + group(Cds-CdsHH) + group(Cds-CdsHH) + radical(Allyl_P)"""),
)

species(
    label = '[CH2]C(C)C(O)=C=C(27849)',
    structure = SMILES('[CH2]C(C)C(O)=C=C'),
    E0 = (82.8285,'kJ/mol'),
    modes = [
        HarmonicOscillator(frequencies=([3615,1277.5,1000,2750,2800,2850,1350,1500,750,1050,1375,1000,350,440,435,1725,540,610,2055,2950,3100,1380,975,1025,1650,3000,3100,440,815,1455,1000,1380,1390,370,380,2900,435,300.753],'cm^-1')),
        HinderedRotor(inertia=(0.248311,'amu*angstrom^2'), symmetry=1, barrier=(15.9729,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.248821,'amu*angstrom^2'), symmetry=1, barrier=(15.9741,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.248272,'amu*angstrom^2'), symmetry=1, barrier=(15.975,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.248867,'amu*angstrom^2'), symmetry=1, barrier=(15.9747,'kJ/mol'), semiclassical=False),
    ],
    spinMultiplicity = 2,
    opticalIsomers = 1,
    molecularWeight = (97.1351,'amu'),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[0.387319,0.0683512,-3.76807e-05,-9.18117e-09,1.15727e-11,10102.3,27.7434], Tmin=(100,'K'), Tmax=(927.654,'K')), NASAPolynomial(coeffs=[18.1921,0.0200571,-5.64059e-06,8.87624e-10,-5.97359e-14,5573.55,-63.4221], Tmin=(927.654,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(82.8285,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(365.837,'J/(mol*K)'), comment="""Thermo group additivity estimation: group(O2s-(Cds-Cd)H) + group(Cs-(Cds-Cdd-Cd)CsCsH) + group(Cs-CsHHH) + group(Cs-CsHHH) + group(Cds-CdsCsOs) + group(Cds-CdsHH) + group(Cdd-CdsCds) + radical(Isobutyl)"""),
)

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
    label = '[CH2]C(=CO)C([CH2])[CH2](14058)',
    structure = SMILES('[CH2]C(=CO)C([CH2])[CH2]'),
    E0 = (269.726,'kJ/mol'),
    modes = [
        HarmonicOscillator(frequencies=([3000,3020,3040,3060,3080,3100,415,440,465,780,815,850,1435,1455,1475,900,1000,1100,3010,987.5,1337.5,450,1655,3615,1277.5,1000,350,440,435,1725,1380,1390,370,380,2900,435,259.757],'cm^-1')),
        HinderedRotor(inertia=(0.0437365,'amu*angstrom^2'), symmetry=1, barrier=(15.8616,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(3.61504,'amu*angstrom^2'), symmetry=1, barrier=(83.1169,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.689977,'amu*angstrom^2'), symmetry=1, barrier=(15.8639,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.043581,'amu*angstrom^2'), symmetry=1, barrier=(15.8593,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.228067,'amu*angstrom^2'), symmetry=1, barrier=(83.139,'kJ/mol'), semiclassical=False),
    ],
    spinMultiplicity = 4,
    opticalIsomers = 1,
    molecularWeight = (97.1351,'amu'),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[-1.38362,0.0868776,-8.43699e-05,4.02653e-08,-7.07372e-12,32661.7,35.1923], Tmin=(100,'K'), Tmax=(1683.29,'K')), NASAPolynomial(coeffs=[21.1258,0.0128709,-1.38163e-07,-3.35871e-10,3.17538e-14,27990.5,-76.49], Tmin=(1683.29,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(269.726,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(361.68,'J/(mol*K)'), comment="""Thermo group additivity estimation: group(O2s-(Cds-Cd)H) + group(Cs-(Cds-Cds)CsCsH) + group(Cs-CsHHH) + group(Cs-CsHHH) + group(Cs-(Cds-Cds)HHH) + group(Cds-CdsCsCs) + longDistanceInteraction_noncyclic(CdCs-ST) + group(Cds-CdsOsH) + radical(Allyl_P) + radical(Isobutyl) + radical(Isobutyl)"""),
)

species(
    label = 'C=[C]C(O)C1CC1(28511)',
    structure = SMILES('C=[C]C(O)C1CC1'),
    E0 = (154.508,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
    molecularWeight = (97.1351,'amu'),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[1.16396,0.0465405,1.54092e-05,-5.48896e-08,2.46562e-11,18699.3,27.7606], Tmin=(100,'K'), Tmax=(974.821,'K')), NASAPolynomial(coeffs=[15.9604,0.0239523,-8.49966e-06,1.58244e-09,-1.15852e-13,14003,-52.5378], Tmin=(974.821,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(154.508,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(369.994,'J/(mol*K)'), comment="""Thermo group additivity estimation: group(O2s-CsH) + group(Cs-CsCsCsH) + group(Cs-CsCsHH) + group(Cs-CsCsHH) + group(Cs-(Cds-Cds)CsOsH) + group(Cds-CdsCsH) + group(Cds-CdsHH) + ring(Cyclopropane) + radical(Cds_S)"""),
)

species(
    label = '[CH2]C1CC(=C)C1O(28512)',
    structure = SMILES('[CH2]C1CC(=C)C1O'),
    E0 = (99.8076,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
    molecularWeight = (97.1351,'amu'),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[10.4846,-0.0115254,0.000128979,-1.24847e-07,2.88446e-11,11655.9,-15.2141], Tmin=(100,'K'), Tmax=(1712.37,'K')), NASAPolynomial(coeffs=[73.9648,0.0285336,-7.10971e-05,1.728e-08,-1.28315e-12,-37697.9,-436.243], Tmin=(1712.37,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(99.8076,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(374.151,'J/(mol*K)'), comment="""Thermo group additivity estimation: group(O2s-CsH) + group(Cs-CsCsCsH) + group(Cs-(Cds-Cds)CsOsH) + group(Cs-(Cds-Cds)CsHH) + group(Cs-CsHHH) + group(Cds-CdsCsCs) + group(Cds-CdsHH) + ring(methylenecyclobutane) + radical(Isobutyl)"""),
)

species(
    label = '[CH2][CH]C(O)[C]=C(27815)',
    structure = SMILES('[CH2][CH]C(O)[C]=C'),
    E0 = (429.881,'kJ/mol'),
    modes = [
        HarmonicOscillator(frequencies=([3615,1277.5,1000,3025,407.5,1350,352.5,2950,3100,1380,975,1025,1650,3000,3100,440,815,1455,1000,1685,370,1380,1390,370,380,2900,435,417.3,417.304],'cm^-1')),
        HinderedRotor(inertia=(0.000968012,'amu*angstrom^2'), symmetry=1, barrier=(0.119627,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.000967935,'amu*angstrom^2'), symmetry=1, barrier=(0.119627,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.000967927,'amu*angstrom^2'), symmetry=1, barrier=(0.119627,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.0703548,'amu*angstrom^2'), symmetry=1, barrier=(8.69422,'kJ/mol'), semiclassical=False),
    ],
    spinMultiplicity = 4,
    opticalIsomers = 1,
    molecularWeight = (83.1085,'amu'),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[1.25612,0.0525465,-4.15588e-05,1.68717e-08,-2.74125e-12,51808.1,30.2053], Tmin=(100,'K'), Tmax=(1469.61,'K')), NASAPolynomial(coeffs=[13.6561,0.018796,-7.1102e-06,1.2445e-09,-8.28526e-14,48163.5,-34.3914], Tmin=(1469.61,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(429.881,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(291.007,'J/(mol*K)'), comment="""Thermo group additivity estimation: group(O2s-CsH) + group(Cs-(Cds-Cds)CsOsH) + group(Cs-CsCsHH) + group(Cs-CsHHH) + group(Cds-CdsCsH) + group(Cds-CdsHH) + radical(CCJCO) + radical(RCCJ) + radical(Cds_S)"""),
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
    label = '[CH2]C([CH2])[CH]O(6286)',
    structure = SMILES('[CH2]C([CH2])[CH]O'),
    E0 = (277.29,'kJ/mol'),
    modes = [
        HarmonicOscillator(frequencies=([3000,3033.33,3066.67,3100,415,465,780,850,1435,1475,900,1100,3025,407.5,1350,352.5,3615,1277.5,1000,1380,1390,370,380,2900,435,180],'cm^-1')),
        HinderedRotor(inertia=(0.00100212,'amu*angstrom^2'), symmetry=1, barrier=(2.97684,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.127566,'amu*angstrom^2'), symmetry=1, barrier=(2.93299,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(2.92283,'amu*angstrom^2'), symmetry=1, barrier=(67.2015,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.128594,'amu*angstrom^2'), symmetry=1, barrier=(2.95662,'kJ/mol'), semiclassical=False),
    ],
    spinMultiplicity = 4,
    opticalIsomers = 1,
    molecularWeight = (71.0978,'amu'),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[1.41279,0.0572235,-6.97008e-05,4.83643e-08,-1.31938e-11,33443,23.7309], Tmin=(100,'K'), Tmax=(1011.97,'K')), NASAPolynomial(coeffs=[9.18421,0.0205716,-6.57753e-06,9.85461e-10,-5.77076e-14,32174,-12.3529], Tmin=(1011.97,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(277.29,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(266.063,'J/(mol*K)'), comment="""Thermo group additivity estimation: group(O2s-CsH) + group(Cs-CsCsCsH) + group(Cs-CsOsHH) + group(Cs-CsHHH) + group(Cs-CsHHH) + radical(Isobutyl) + radical(Isobutyl) + radical(CCsJOH)"""),
)

species(
    label = '[CH]C([CH2])C(O)[C]=C(28513)',
    structure = SMILES('[CH]C([CH2])C(O)[C]=C'),
    E0 = (645.209,'kJ/mol'),
    modes = [
        HarmonicOscillator(frequencies=([3000,3100,440,815,1455,1000,2950,3100,1380,975,1025,1650,3615,1277.5,1000,1685,370,1380,1383.33,1386.67,1390,370,373.333,376.667,380,2800,3000,430,440,200,800,1066.67,1333.33,1600],'cm^-1')),
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
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[0.508411,0.0724031,-7.06242e-05,3.67398e-08,-7.64621e-12,77730.1,33.3458], Tmin=(100,'K'), Tmax=(1165.15,'K')), NASAPolynomial(coeffs=[14.4029,0.0247025,-9.21455e-06,1.60265e-09,-1.06979e-13,74492.3,-35.8108], Tmin=(1165.15,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(645.209,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(336.736,'J/(mol*K)'), comment="""Thermo group additivity estimation: group(O2s-CsH) + group(Cs-CsCsCsH) + group(Cs-(Cds-Cds)CsOsH) + group(Cs-CsHHH) + group(Cs-CsHHH) + group(Cds-CdsCsH) + group(Cds-CdsHH) + radical(CCJ2_triplet) + radical(Cds_S) + radical(Isobutyl)"""),
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
    E0 = (402.076,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS2',
    E0 = (484.296,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS3',
    E0 = (512.14,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS4',
    E0 = (555.139,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS5',
    E0 = (530.128,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS6',
    E0 = (402.076,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS7',
    E0 = (480.377,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS8',
    E0 = (530.864,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS9',
    E0 = (543.564,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS10',
    E0 = (584.657,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS11',
    E0 = (516.767,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS12',
    E0 = (520.015,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS13',
    E0 = (544.053,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS14',
    E0 = (544.053,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS15',
    E0 = (485.756,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS16',
    E0 = (477.13,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS17',
    E0 = (743.705,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS18',
    E0 = (851.336,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS19',
    E0 = (670.965,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS20',
    E0 = (766.455,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS21',
    E0 = (712.452,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS22',
    E0 = (860.977,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS23',
    E0 = (424.937,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS24',
    E0 = (491.045,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS25',
    E0 = (480.323,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS26',
    E0 = (559.394,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS27',
    E0 = (496.551,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS28',
    E0 = (410.165,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS29',
    E0 = (410.36,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS30',
    E0 = (845.567,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS31',
    E0 = (911.857,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS32',
    E0 = (857.013,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

reaction(
    label = 'reaction1',
    reactants = ['[CH2]C([CH2])C(O)[C]=C(14106)'],
    products = ['[CH2]C=C(87)', 'C=C=CO(12571)'],
    transitionState = 'TS1',
    kinetics = Arrhenius(A=(5e+12,'s^-1'), n=0, Ea=(0,'kJ/mol'), T0=(1,'K'), Tmin=(300,'K'), Tmax=(1500,'K'), comment="""Exact match found for rate rule [RJJ]
Euclidian distance = 0
family: 1,4_Linear_birad_scission"""),
)

reaction(
    label = 'reaction2',
    reactants = ['H(8)', '[CH2]C(=C)C(O)[C]=C(18301)'],
    products = ['[CH2]C([CH2])C(O)[C]=C(14106)'],
    transitionState = 'TS2',
    kinetics = Arrhenius(A=(170.395,'m^3/(mol*s)'), n=1.5621, Ea=(11.2886,'kJ/mol'), T0=(1,'K'), comment="""Estimated using an average for rate rule [Cds_Cds;HJ]
Euclidian distance = 0
family: R_Addition_MultipleBond"""),
)

reaction(
    label = 'reaction3',
    reactants = ['H(8)', '[CH2]C([CH2])C(O)=C=C(28507)'],
    products = ['[CH2]C([CH2])C(O)[C]=C(14106)'],
    transitionState = 'TS3',
    kinetics = Arrhenius(A=(169.619,'m^3/(mol*s)'), n=1.605, Ea=(12.4249,'kJ/mol'), T0=(1,'K'), comment="""Estimated using an average for rate rule [Cds_Ca;HJ]
Euclidian distance = 0
family: R_Addition_MultipleBond"""),
)

reaction(
    label = 'reaction4',
    reactants = ['H(8)', 'C#CC(O)C([CH2])[CH2](22792)'],
    products = ['[CH2]C([CH2])C(O)[C]=C(14106)'],
    transitionState = 'TS4',
    kinetics = Arrhenius(A=(1.255e+11,'cm^3/(mol*s)'), n=1.005, Ea=(13.1503,'kJ/mol'), T0=(1,'K'), comment="""From training reaction 138 used for Ct-H_Ct-Cs;HJ
Exact match found for rate rule [Ct-H_Ct-Cs;HJ]
Euclidian distance = 0
family: R_Addition_MultipleBond"""),
)

reaction(
    label = 'reaction5',
    reactants = ['CH2(T)(28)', 'C=[C]C(O)C=C(27708)'],
    products = ['[CH2]C([CH2])C(O)[C]=C(14106)'],
    transitionState = 'TS5',
    kinetics = Arrhenius(A=(4.17e+07,'cm^3/(mol*s)'), n=1.64, Ea=(0,'kJ/mol'), T0=(1,'K'), Tmin=(300,'K'), Tmax=(1500,'K'), comment="""Estimated using template [Cds-CsH_Cds-HH;Y_1centerbirad] for rate rule [Cds-Cs\O2s/H_Cds-HH;CH2_triplet]
Euclidian distance = 1.41421356237
family: R_Addition_MultipleBond
Ea raised from -5.9 to 0 kJ/mol."""),
)

reaction(
    label = 'reaction6',
    reactants = ['[CH2]C=C(87)', '[CH2][C]=CO(18753)'],
    products = ['[CH2]C([CH2])C(O)[C]=C(14106)'],
    transitionState = 'TS6',
    kinetics = Arrhenius(A=(0.00179511,'m^3/(mol*s)'), n=2.50446, Ea=(57.781,'kJ/mol'), T0=(1,'K'), comment="""Estimated using an average for rate rule [Cds_Cds;CJ]
Euclidian distance = 0
family: R_Addition_MultipleBond
Ea raised from 53.4 to 57.8 kJ/mol to match endothermicity of reaction."""),
)

reaction(
    label = 'reaction7',
    reactants = ['[CH2][CH][CH2](6136)', 'C=C=CO(12571)'],
    products = ['[CH2]C([CH2])C(O)[C]=C(14106)'],
    transitionState = 'TS7',
    kinetics = Arrhenius(A=(0.00401797,'m^3/(mol*s)'), n=2.41733, Ea=(22.1495,'kJ/mol'), T0=(1,'K'), comment="""Estimated using an average for rate rule [Cds_Ca;CJ]
Euclidian distance = 0
family: R_Addition_MultipleBond"""),
)

reaction(
    label = 'reaction8',
    reactants = ['OH(D)(132)', '[CH2]C([CH2])C=C=C(17716)'],
    products = ['[CH2]C([CH2])C(O)[C]=C(14106)'],
    transitionState = 'TS8',
    kinetics = Arrhenius(A=(986500,'cm^3/(mol*s)'), n=2.037, Ea=(0,'kJ/mol'), T0=(1,'K'), comment="""Estimated using template [Cds_Ca;OJ_pri] for rate rule [Cds-CsH_Ca;OJ_pri]
Euclidian distance = 1.0
family: R_Addition_MultipleBond
Ea raised from -6.0 to 0 kJ/mol."""),
)

reaction(
    label = 'reaction9',
    reactants = ['[CH2]C([CH2])C(O)[C]=C(14106)'],
    products = ['[CH2][C](C)C(O)[C]=C(27854)'],
    transitionState = 'TS9',
    kinetics = Arrhenius(A=(34962.3,'s^-1'), n=2.56136, Ea=(141.488,'kJ/mol'), T0=(1,'K'), comment="""Estimated using an average for rate rule [R2H_S;C_rad_out_2H;XH_out]
Euclidian distance = 0
Multiplied by reaction path degeneracy 2.0
family: intra_H_migration"""),
)

reaction(
    label = 'reaction10',
    reactants = ['[CH2]C([CH2])C(O)[C]=C(14106)'],
    products = ['[CH2]C=C(O)C([CH2])[CH2](14104)'],
    transitionState = 'TS10',
    kinetics = Arrhenius(A=(3.677e+10,'s^-1'), n=0.839, Ea=(182.581,'kJ/mol'), T0=(1,'K'), comment="""Estimated using template [R2H_S;Cd_rad_out_Cd;Cs_H_out_noH] for rate rule [R2H_S;Cd_rad_out_Cd;Cs_H_out_NDMustO]
Euclidian distance = 2.0
family: intra_H_migration"""),
)

reaction(
    label = 'reaction11',
    reactants = ['[CH]=CC(O)C([CH2])[CH2](14107)'],
    products = ['[CH2]C([CH2])C(O)[C]=C(14106)'],
    transitionState = 'TS11',
    kinetics = Arrhenius(A=(1.08e+06,'s^-1'), n=1.99, Ea=(105.437,'kJ/mol'), T0=(1,'K'), comment="""From training reaction 17 used for R2H_D;Cd_rad_out_singleH;Cd_H_out_singleNd
Exact match found for rate rule [R2H_D;Cd_rad_out_singleH;Cd_H_out_singleNd]
Euclidian distance = 0
family: intra_H_migration"""),
)

reaction(
    label = 'reaction12',
    reactants = ['[CH2]C([CH2])C(O)[C]=C(14106)'],
    products = ['[CH2][C]=C(O)C([CH2])C(27855)'],
    transitionState = 'TS12',
    kinetics = Arrhenius(A=(166690,'s^-1'), n=2.17519, Ea=(117.939,'kJ/mol'), T0=(1,'K'), comment="""Estimated using an average for rate rule [R3H_SS_Cs;C_rad_out_2H;XH_out]
Euclidian distance = 0
Multiplied by reaction path degeneracy 2.0
family: intra_H_migration"""),
)

reaction(
    label = 'reaction14',
    reactants = ['[CH2]C([CH2])C(O)[C]=C(14106)'],
    products = ['[CH2]C([CH2])C([O])C=C(12805)'],
    transitionState = 'TS13',
    kinetics = Arrhenius(A=(2.4115e+09,'s^-1'), n=1.00333, Ea=(141.977,'kJ/mol'), T0=(1,'K'), comment="""Estimated using template [R3H_SS_Cs;Cd_rad_out_Cd;XH_out] for rate rule [R3H_SS_Cs;Cd_rad_out_Cd;O_H_out]
Euclidian distance = 1.0
family: intra_H_migration"""),
)

reaction(
    label = 'reaction14',
    reactants = ['[CH2]C([CH2])C(O)[C]=C(14106)'],
    products = ['[CH2][C]([CH2])C(O)C=C(14105)'],
    transitionState = 'TS14',
    kinetics = Arrhenius(A=(2.4115e+09,'s^-1'), n=1.00333, Ea=(141.977,'kJ/mol'), T0=(1,'K'), comment="""Estimated using an average for rate rule [R3H_SS_Cs;Cd_rad_out_Cd;XH_out]
Euclidian distance = 0
family: intra_H_migration"""),
)

reaction(
    label = 'reaction15',
    reactants = ['[CH2]C([CH2])C(O)[C]=C(14106)'],
    products = ['[CH2]C(C)C([O])[C]=C(13089)'],
    transitionState = 'TS15',
    kinetics = Arrhenius(A=(1.72e-08,'s^-1'), n=5.55, Ea=(83.68,'kJ/mol'), T0=(1,'K'), Tmin=(300,'K'), Tmax=(2500,'K'), comment="""From training reaction 340 used for R4H_SSS;C_rad_out_2H;O_H_out
Exact match found for rate rule [R4H_SSS;C_rad_out_2H;O_H_out]
Euclidian distance = 0
Multiplied by reaction path degeneracy 2.0
family: intra_H_migration"""),
)

reaction(
    label = 'reaction16',
    reactants = ['[CH]=[C]C(O)C([CH2])C(27856)'],
    products = ['[CH2]C([CH2])C(O)[C]=C(14106)'],
    transitionState = 'TS16',
    kinetics = Arrhenius(A=(408000,'s^-1'), n=1.9199, Ea=(33.0402,'kJ/mol'), T0=(1,'K'), Tmin=(300,'K'), Tmax=(1600,'K'), comment="""Estimated using template [R5Hall;Cd_rad_out_singleH;Cs_H_out_2H] for rate rule [R5HJ_1;Cd_rad_out_singleH;Cs_H_out_2H]
Euclidian distance = 1.0
Multiplied by reaction path degeneracy 3.0
family: intra_H_migration"""),
)

reaction(
    label = 'reaction17',
    reactants = ['OH(D)(132)', '[CH2][C]=CC([CH2])[CH2](17719)'],
    products = ['[CH2]C([CH2])C(O)[C]=C(14106)'],
    transitionState = 'TS17',
    kinetics = Arrhenius(A=(3.05166e+07,'m^3/(mol*s)'), n=0.045, Ea=(0.1046,'kJ/mol'), T0=(1,'K'), comment="""Estimated using an average for rate rule [O_pri_rad;Y_rad]
Euclidian distance = 0
family: R_Recombination"""),
)

reaction(
    label = 'reaction18',
    reactants = ['H(8)', '[CH2]C([CH2])C([O])[C]=C(14110)'],
    products = ['[CH2]C([CH2])C(O)[C]=C(14106)'],
    transitionState = 'TS18',
    kinetics = Arrhenius(A=(5.00518e+06,'m^3/(mol*s)'), n=0.282325, Ea=(7.09479,'kJ/mol'), T0=(1,'K'), comment="""Estimated using average of templates [Y_rad;O_rad/NonDe] + [H_rad;O_sec_rad] for rate rule [H_rad;O_rad/NonDe]
Euclidian distance = 1.0
family: R_Recombination"""),
)

reaction(
    label = 'reaction19',
    reactants = ['[CH2][CH][CH2](6136)', '[CH2][C]=CO(18753)'],
    products = ['[CH2]C([CH2])C(O)[C]=C(14106)'],
    transitionState = 'TS19',
    kinetics = Arrhenius(A=(1.9789e+07,'m^3/(mol*s)'), n=-0.126319, Ea=(0,'kJ/mol'), T0=(1,'K'), comment="""Estimated using an average for rate rule [Y_rad;Y_rad]
Euclidian distance = 0
family: R_Recombination
Ea raised from -15.6 to -15.6 kJ/mol.
Ea raised from -15.6 to 0 kJ/mol."""),
)

reaction(
    label = 'reaction20',
    reactants = ['H(8)', '[CH2][C]([CH2])C(O)[C]=C(28508)'],
    products = ['[CH2]C([CH2])C(O)[C]=C(14106)'],
    transitionState = 'TS20',
    kinetics = Arrhenius(A=(4.34078e+06,'m^3/(mol*s)'), n=0.278577, Ea=(0,'kJ/mol'), T0=(1,'K'), comment="""Estimated using an average for rate rule [Y_rad;H_rad]
Euclidian distance = 0
family: R_Recombination
Ea raised from -1.4 to 0 kJ/mol."""),
)

reaction(
    label = 'reaction21',
    reactants = ['H(8)', '[CH2][C]=C(O)C([CH2])[CH2](28509)'],
    products = ['[CH2]C([CH2])C(O)[C]=C(14106)'],
    transitionState = 'TS21',
    kinetics = Arrhenius(A=(4.34078e+06,'m^3/(mol*s)'), n=0.278577, Ea=(0,'kJ/mol'), T0=(1,'K'), comment="""Estimated using an average for rate rule [Y_rad;H_rad]
Euclidian distance = 0
family: R_Recombination
Ea raised from -1.4 to 0 kJ/mol."""),
)

reaction(
    label = 'reaction22',
    reactants = ['H(8)', '[CH]=[C]C(O)C([CH2])[CH2](28510)'],
    products = ['[CH2]C([CH2])C(O)[C]=C(14106)'],
    transitionState = 'TS22',
    kinetics = Arrhenius(A=(4.34078e+06,'m^3/(mol*s)'), n=0.278577, Ea=(0,'kJ/mol'), T0=(1,'K'), comment="""Estimated using an average for rate rule [Y_rad;H_rad]
Euclidian distance = 0
family: R_Recombination
Ea raised from -1.4 to 0 kJ/mol."""),
)

reaction(
    label = 'reaction23',
    reactants = ['[CH2]C([CH2])C(O)[C]=C(14106)'],
    products = ['C=[C]C(O)C(=C)C(27848)'],
    transitionState = 'TS23',
    kinetics = Arrhenius(A=(3.898e+11,'s^-1'), n=0.486, Ea=(22.8614,'kJ/mol'), T0=(1,'K'), comment="""Estimated using an average for rate rule [R2radExo;Y_rad;XH_Rrad]
Euclidian distance = 0
Multiplied by reaction path degeneracy 2.0
family: Intra_Disproportionation"""),
)

reaction(
    label = 'reaction24',
    reactants = ['[CH2]C([CH2])C(O)[C]=C(14106)'],
    products = ['[CH2]C(=C)C(O)C=C(14114)'],
    transitionState = 'TS24',
    kinetics = Arrhenius(A=(2.6374e+09,'s^-1'), n=0.37, Ea=(88.9686,'kJ/mol'), T0=(1,'K'), comment="""Estimated using average of templates [R3;Y_rad_De;XH_Rrad] + [R3radExo;Y_rad;XH_Rrad] for rate rule [R3radExo;Y_rad_De;XH_Rrad]
Euclidian distance = 1.0
family: Intra_Disproportionation"""),
)

reaction(
    label = 'reaction25',
    reactants = ['[CH2]C([CH2])C(O)[C]=C(14106)'],
    products = ['[CH2]C(C)C(O)=C=C(27849)'],
    transitionState = 'TS25',
    kinetics = Arrhenius(A=(4.00798e+09,'s^-1'), n=0.37, Ea=(78.2471,'kJ/mol'), T0=(1,'K'), comment="""Estimated using average of templates [R3;Y_rad;XH_Rrad_De] + [R3radExo;Y_rad;XH_Rrad] for rate rule [R3radExo;Y_rad;XH_Rrad_De]
Euclidian distance = 1.0
Multiplied by reaction path degeneracy 2.0
family: Intra_Disproportionation"""),
)

reaction(
    label = 'reaction26',
    reactants = ['[CH2]C([CH2])C(O)[C]=C(14106)'],
    products = ['[CH2][CH]CC(O)[C]=C(14200)'],
    transitionState = 'TS26',
    kinetics = Arrhenius(A=(2.66e+08,'s^-1'), n=1.36, Ea=(157.318,'kJ/mol'), T0=(1,'K'), comment="""Estimated using an average for rate rule [cCsCJ;CsJ-HH;C]
Euclidian distance = 0
Multiplied by reaction path degeneracy 2.0
family: 1,2_shiftC"""),
)

reaction(
    label = 'reaction27',
    reactants = ['[CH2]C([CH2])C(O)[C]=C(14106)'],
    products = ['[CH2]C(=CO)C([CH2])[CH2](14058)'],
    transitionState = 'TS27',
    kinetics = Arrhenius(A=(8.66e+11,'s^-1'), n=0.438, Ea=(94.4747,'kJ/mol'), T0=(1,'K'), comment="""Estimated using an average for rate rule [cCs(-HR!H)CJ;CdsJ;C]
Euclidian distance = 0
family: 1,2_shiftC"""),
)

reaction(
    label = 'reaction28',
    reactants = ['[CH2]C([CH2])C(O)[C]=C(14106)'],
    products = ['C=[C]C(O)C1CC1(28511)'],
    transitionState = 'TS28',
    kinetics = Arrhenius(A=(7.38971e+10,'s^-1'), n=0.0476667, Ea=(8.08907,'kJ/mol'), T0=(1,'K'), comment="""Estimated using template [Rn;C_rad_out_2H;Cpri_rad_out_2H] for rate rule [R3_SS;C_rad_out_2H;Cpri_rad_out_2H]
Euclidian distance = 2.0
family: Birad_recombination"""),
)

reaction(
    label = 'reaction29',
    reactants = ['[CH2]C([CH2])C(O)[C]=C(14106)'],
    products = ['[CH2]C1CC(=C)C1O(28512)'],
    transitionState = 'TS29',
    kinetics = Arrhenius(A=(3.24e+12,'s^-1'), n=-0.305, Ea=(8.28432,'kJ/mol'), T0=(1,'K'), Tmin=(600,'K'), Tmax=(2000,'K'), comment="""Estimated using an average for rate rule [R4_SSS;C_rad_out_2H;Ypri_rad_out]
Euclidian distance = 0
Multiplied by reaction path degeneracy 2.0
family: Birad_recombination"""),
)

reaction(
    label = 'reaction30',
    reactants = ['CH2(T)(28)', '[CH2][CH]C(O)[C]=C(27815)'],
    products = ['[CH2]C([CH2])C(O)[C]=C(14106)'],
    transitionState = 'TS30',
    kinetics = Arrhenius(A=(1.14854e+06,'m^3/(mol*s)'), n=0.575199, Ea=(34.3157,'kJ/mol'), T0=(1,'K'), comment="""Estimated using an average for rate rule [Y_rad;Birad]
Euclidian distance = 0
family: Birad_R_Recombination"""),
)

reaction(
    label = 'reaction31',
    reactants = ['[C]=C(584)', '[CH2]C([CH2])[CH]O(6286)'],
    products = ['[CH2]C([CH2])C(O)[C]=C(14106)'],
    transitionState = 'TS31',
    kinetics = Arrhenius(A=(1.14854e+06,'m^3/(mol*s)'), n=0.575199, Ea=(34.3157,'kJ/mol'), T0=(1,'K'), comment="""Estimated using template [Y_rad;Birad] for rate rule [C_rad/H/CsO;Birad]
Euclidian distance = 4.0
family: Birad_R_Recombination"""),
)

reaction(
    label = 'reaction32',
    reactants = ['H(8)', '[CH]C([CH2])C(O)[C]=C(28513)'],
    products = ['[CH2]C([CH2])C(O)[C]=C(14106)'],
    transitionState = 'TS32',
    kinetics = Arrhenius(A=(1e+07,'m^3/(mol*s)'), n=0, Ea=(0,'kJ/mol'), T0=(1,'K'), comment="""Estimated using an average for rate rule [H_rad;Birad]
Euclidian distance = 0
family: Birad_R_Recombination"""),
)

network(
    label = '5140',
    isomers = [
        '[CH2]C([CH2])C(O)[C]=C(14106)',
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
    label = '5140',
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

