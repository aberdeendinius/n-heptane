species(
    label = '[O]OC([O])C[C]=CO(15241)',
    structure = SMILES('[O]OC([O])C[C]=CO'),
    E0 = (101.684,'kJ/mol'),
    modes = [
        HarmonicOscillator(frequencies=([492.5,1135,1000,3615,1277.5,1000,2750,2850,1437.5,1250,1305,750,350,1380,1390,370,380,2900,435,1685,370,3010,987.5,1337.5,450,1655,297.553,297.553,297.553],'cm^-1')),
        HinderedRotor(inertia=(0.203209,'amu*angstrom^2'), symmetry=1, barrier=(12.7673,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.203209,'amu*angstrom^2'), symmetry=1, barrier=(12.7673,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.203209,'amu*angstrom^2'), symmetry=1, barrier=(12.7673,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.203209,'amu*angstrom^2'), symmetry=1, barrier=(12.7673,'kJ/mol'), semiclassical=False),
    ],
    spinMultiplicity = 4,
    opticalIsomers = 1,
    molecularWeight = (117.08,'amu'),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[0.336584,0.0833465,-0.000109578,7.44152e-08,-1.99569e-11,12359.3,31.8131], Tmin=(100,'K'), Tmax=(914.415,'K')), NASAPolynomial(coeffs=[14.1969,0.0227154,-1.01185e-05,1.90212e-09,-1.3178e-13,9824.45,-33.8148], Tmin=(914.415,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(101.684,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(291.007,'J/(mol*K)'), comment="""Thermo group additivity estimation: group(O2s-OsCs) + group(O2s-CsH) + group(O2s-(Cds-Cd)H) + group(O2s-OsH) + group(Cs-(Cds-Cds)CsHH) + group(Cs-CsOsOsH) + group(Cds-CdsCsH) + group(Cds-CdsOsH) + radical(CCOJ) + radical(Cds_S) + radical(ROOJ)"""),
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
    label = '[O]OC([O])C=C=CO(28885)',
    structure = SMILES('[O]OC([O])C=C=CO'),
    E0 = (27.0661,'kJ/mol'),
    modes = [
        HarmonicOscillator(frequencies=([492.5,1135,1000,540,610,2055,2995,3025,975,1000,1300,1375,400,500,1630,1680,3615,1277.5,1000,1380,1390,370,380,2900,435,180,180],'cm^-1')),
        HinderedRotor(inertia=(0.825598,'amu*angstrom^2'), symmetry=1, barrier=(18.9821,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.82527,'amu*angstrom^2'), symmetry=1, barrier=(18.9746,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.825322,'amu*angstrom^2'), symmetry=1, barrier=(18.9758,'kJ/mol'), semiclassical=False),
    ],
    spinMultiplicity = 3,
    opticalIsomers = 1,
    molecularWeight = (116.072,'amu'),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[0.397019,0.0822615,-0.000112556,7.83299e-08,-2.1408e-11,3382.41,29.0697], Tmin=(100,'K'), Tmax=(899.453,'K')), NASAPolynomial(coeffs=[14.3277,0.0203078,-9.23372e-06,1.74622e-09,-1.21163e-13,876.485,-36.6612], Tmin=(899.453,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(27.0661,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(270.22,'J/(mol*K)'), comment="""Thermo group additivity estimation: group(O2s-OsCs) + group(O2s-CsH) + group(O2s-(Cds-Cd)H) + group(O2s-OsH) + group(Cs-CsOsOsH) + group(Cds-CdsCsH) + group(Cds-CdsOsH) + group(Cdd-CdsCds) + radical(ROOJ) + radical(CCOJ)"""),
)

species(
    label = '[O]OC(=O)C[C]=CO(28886)',
    structure = SMILES('[O]OC(=O)C[C]=CO'),
    E0 = (-137.548,'kJ/mol'),
    modes = [
        HarmonicOscillator(frequencies=([1685,370,492.5,1135,1000,3010,987.5,1337.5,450,1655,3615,1277.5,1000,2750,2850,1437.5,1250,1305,750,350,200,800,1000,1200,1400,1600],'cm^-1')),
        HinderedRotor(inertia=(0.156089,'amu*angstrom^2'), symmetry=1, barrier=(3.5888,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.156089,'amu*angstrom^2'), symmetry=1, barrier=(3.5888,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.156089,'amu*angstrom^2'), symmetry=1, barrier=(3.5888,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.156089,'amu*angstrom^2'), symmetry=1, barrier=(3.5888,'kJ/mol'), semiclassical=False),
    ],
    spinMultiplicity = 3,
    opticalIsomers = 1,
    molecularWeight = (116.072,'amu'),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[0.593112,0.056479,-1.07975e-05,-4.12229e-08,2.32232e-11,-16403.9,28.5231], Tmin=(100,'K'), Tmax=(963.333,'K')), NASAPolynomial(coeffs=[23.5107,0.00590273,-1.46608e-06,3.61326e-10,-3.62135e-14,-22888,-91.9221], Tmin=(963.333,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(-137.548,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(266.063,'J/(mol*K)'), comment="""Thermo group additivity estimation: group(O2s-O2s(Cds-O2d)) + group(O2s-(Cds-Cd)H) + group(O2s-OsH) + group(Cs-(Cds-O2d)(Cds-Cds)HH) + group(Cds-CdsCsH) + group(Cds-OdCsOs) + group(Cds-CdsOsH) + radical(C(=O)OOJ) + radical(Cds_S)"""),
)

species(
    label = '[O]OC([O])CC#CO(28887)',
    structure = SMILES('[O]OC([O])CC#CO'),
    E0 = (97.2722,'kJ/mol'),
    modes = [
        HarmonicOscillator(frequencies=([492.5,1135,1000,1380,1390,370,380,2900,435,3615,1277.5,1000,2100,2250,500,550,2750,2850,1437.5,1250,1305,750,350,180,180,1215.5],'cm^-1')),
        HinderedRotor(inertia=(0.201893,'amu*angstrom^2'), symmetry=1, barrier=(4.64191,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.204225,'amu*angstrom^2'), symmetry=1, barrier=(4.69553,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.20391,'amu*angstrom^2'), symmetry=1, barrier=(4.68829,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(2.46811,'amu*angstrom^2'), symmetry=1, barrier=(56.7468,'kJ/mol'), semiclassical=False),
    ],
    spinMultiplicity = 3,
    opticalIsomers = 1,
    molecularWeight = (116.072,'amu'),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[0.905258,0.0770762,-0.000126582,1.15713e-07,-4.10491e-11,11801.9,30.6094], Tmin=(100,'K'), Tmax=(847.193,'K')), NASAPolynomial(coeffs=[6.28528,0.0331922,-1.61594e-05,3.06916e-09,-2.09807e-13,11553.6,9.46054], Tmin=(847.193,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(97.2722,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(266.063,'J/(mol*K)'), comment="""Thermo group additivity estimation: group(O2s-OsCs) + group(O2s-CsH) + group(O2s-OsH) + group(O2s-CtH) + group(Cs-CsOsOsH) + group(Cs-CtCsHH) + group(Ct-CtCs) + group(Ct-CtOs) + radical(CCOJ) + radical(ROOJ)"""),
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
    label = 'C#CCC([O])O[O](18040)',
    structure = SMILES('C#CCC([O])O[O]'),
    E0 = (238.738,'kJ/mol'),
    modes = [
        HarmonicOscillator(frequencies=([492.5,1135,1000,750,770,3400,2100,1380,1390,370,380,2900,435,2175,525,2750,2850,1437.5,1250,1305,750,350,180,180],'cm^-1')),
        HinderedRotor(inertia=(0.0049375,'amu*angstrom^2'), symmetry=1, barrier=(7.97069,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.346169,'amu*angstrom^2'), symmetry=1, barrier=(7.95911,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(2.60626,'amu*angstrom^2'), symmetry=1, barrier=(59.9231,'kJ/mol'), semiclassical=False),
    ],
    spinMultiplicity = 3,
    opticalIsomers = 1,
    molecularWeight = (100.073,'amu'),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[1.26195,0.0680287,-0.000110201,9.93603e-08,-3.45787e-11,28804.7,25.9926], Tmin=(100,'K'), Tmax=(867.362,'K')), NASAPolynomial(coeffs=[6.03973,0.0290494,-1.34851e-05,2.49861e-09,-1.6799e-13,28613.3,7.29692], Tmin=(867.362,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(238.738,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(245.277,'J/(mol*K)'), comment="""Thermo group additivity estimation: group(O2s-OsCs) + group(O2s-CsH) + group(O2s-OsH) + group(Cs-CsOsOsH) + group(Cs-CtCsHH) + group(Ct-CtCs) + group(Ct-CtH) + radical(ROOJ) + radical(CCOJ)"""),
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
    label = 'O=CC[C]=CO(27745)',
    structure = SMILES('O=CC[C]=CO'),
    E0 = (-64.0896,'kJ/mol'),
    modes = [
        HarmonicOscillator(frequencies=([1685,370,2750,2850,1437.5,1250,1305,750,350,3010,987.5,1337.5,450,1655,3615,1277.5,1000,2782.5,750,1395,475,1775,1000,180],'cm^-1')),
        HinderedRotor(inertia=(0.935575,'amu*angstrom^2'), symmetry=1, barrier=(21.5107,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.93612,'amu*angstrom^2'), symmetry=1, barrier=(21.5232,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.935243,'amu*angstrom^2'), symmetry=1, barrier=(21.5031,'kJ/mol'), semiclassical=False),
    ],
    spinMultiplicity = 2,
    opticalIsomers = 1,
    molecularWeight = (85.0813,'amu'),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[1.49226,0.0417242,-4.4337e-06,-3.09691e-08,1.61113e-11,-7605.93,20.9994], Tmin=(100,'K'), Tmax=(988.339,'K')), NASAPolynomial(coeffs=[17.1066,0.00994677,-3.8861e-06,8.23847e-10,-6.62225e-14,-12226.8,-61.91], Tmin=(988.339,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(-64.0896,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(245.277,'J/(mol*K)'), comment="""Thermo group additivity estimation: group(O2s-(Cds-Cd)H) + group(Cs-(Cds-O2d)(Cds-Cds)HH) + group(Cds-CdsCsH) + group(Cds-CdsOsH) + group(Cds-OdCsH) + radical(Cds_S)"""),
)

species(
    label = '[O]OC([O])[CH]C=CO(15244)',
    structure = SMILES('[O]OC([O])[CH]C=CO'),
    E0 = (-19.2414,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
    molecularWeight = (117.08,'amu'),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[0.0420311,0.0840781,-9.77112e-05,5.60877e-08,-1.25349e-11,-2169.06,29.5891], Tmin=(100,'K'), Tmax=(1099.16,'K')), NASAPolynomial(coeffs=[18.1946,0.0180179,-7.55978e-06,1.40843e-09,-9.82826e-14,-6159.56,-59.7029], Tmin=(1099.16,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(-19.2414,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(291.007,'J/(mol*K)'), comment="""Thermo group additivity estimation: group(O2s-OsCs) + group(O2s-CsH) + group(O2s-(Cds-Cd)H) + group(O2s-OsH) + group(Cs-(Cds-Cds)CsHH) + group(Cs-CsOsOsH) + group(Cds-CdsCsH) + group(Cds-CdsOsH) + radical(C=CCJCO) + radical(ROOJ) + radical(CCOJ)"""),
)

species(
    label = '[O]O[C](O)C[C]=CO(28888)',
    structure = SMILES('[O]O[C](O)C[C]=CO'),
    E0 = (81.2249,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
    molecularWeight = (117.08,'amu'),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[-0.172773,0.0891379,-0.000116468,7.38512e-08,-1.7951e-11,9921.71,33.4593], Tmin=(100,'K'), Tmax=(1020.69,'K')), NASAPolynomial(coeffs=[19.2568,0.0129956,-4.57039e-06,7.65295e-10,-5.0041e-14,5955.38,-60.675], Tmin=(1020.69,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(81.2249,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(286.849,'J/(mol*K)'), comment="""Thermo group additivity estimation: group(O2s-OsCs) + group(O2s-CsH) + group(O2s-(Cds-Cd)H) + group(O2s-OsH) + group(Cs-(Cds-Cds)CsHH) + group(Cs-CsOsOsH) + group(Cds-CdsCsH) + group(Cds-CdsOsH) + radical(ROOJ) + radical(Cds_S) + radical(Cs_P)"""),
)

species(
    label = '[O]OC([O])CC=[C]O(15239)',
    structure = SMILES('[O]OC([O])CC=[C]O'),
    E0 = (103.586,'kJ/mol'),
    modes = [
        HarmonicOscillator(frequencies=([492.5,1135,1000,3615,1277.5,1000,2750,2850,1437.5,1250,1305,750,350,1380,1390,370,380,2900,435,1685,370,3010,987.5,1337.5,450,1655,256.872,256.874,2177.62],'cm^-1')),
        HinderedRotor(inertia=(0.272539,'amu*angstrom^2'), symmetry=1, barrier=(12.7585,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.272413,'amu*angstrom^2'), symmetry=1, barrier=(12.7586,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.272535,'amu*angstrom^2'), symmetry=1, barrier=(12.7585,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.272507,'amu*angstrom^2'), symmetry=1, barrier=(12.7586,'kJ/mol'), semiclassical=False),
    ],
    spinMultiplicity = 4,
    opticalIsomers = 1,
    molecularWeight = (117.08,'amu'),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[0.534555,0.0810615,-0.000113781,8.69829e-08,-2.67046e-11,12578.8,33.7557], Tmin=(100,'K'), Tmax=(796.796,'K')), NASAPolynomial(coeffs=[11.2266,0.0273804,-1.27131e-05,2.41134e-09,-1.66751e-13,10875.1,-15.3972], Tmin=(796.796,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(103.586,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(291.007,'J/(mol*K)'), comment="""Thermo group additivity estimation: group(O2s-OsCs) + group(O2s-CsH) + group(O2s-(Cds-Cd)H) + group(O2s-OsH) + group(Cs-(Cds-Cds)CsHH) + group(Cs-CsOsOsH) + group(Cds-CdsCsH) + group(Cds-CdsOsH) + radical(CCOJ) + radical(C=CJO) + radical(ROOJ)"""),
)

species(
    label = '[O][C](C[C]=CO)OO(28889)',
    structure = SMILES('[O][C](C[C]=CO)OO'),
    E0 = (154.925,'kJ/mol'),
    modes = [
        HarmonicOscillator(frequencies=([360,370,350,3615,1277.5,1000,2750,2850,1437.5,1250,1305,750,350,3615,1310,387.5,850,1000,1685,370,3010,987.5,1337.5,450,1655,200,800,1600],'cm^-1')),
        HinderedRotor(inertia=(0.156089,'amu*angstrom^2'), symmetry=1, barrier=(3.5888,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.156089,'amu*angstrom^2'), symmetry=1, barrier=(3.5888,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.156089,'amu*angstrom^2'), symmetry=1, barrier=(3.5888,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.156089,'amu*angstrom^2'), symmetry=1, barrier=(3.5888,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.156089,'amu*angstrom^2'), symmetry=1, barrier=(3.5888,'kJ/mol'), semiclassical=False),
    ],
    spinMultiplicity = 4,
    opticalIsomers = 1,
    molecularWeight = (117.08,'amu'),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[0.167253,0.0869674,-0.000112724,7.3778e-08,-1.90092e-11,18769,32.3518], Tmin=(100,'K'), Tmax=(951.131,'K')), NASAPolynomial(coeffs=[15.6093,0.0220257,-1.03068e-05,1.99171e-09,-1.40564e-13,15831.5,-41.3735], Tmin=(951.131,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(154.925,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(286.849,'J/(mol*K)'), comment="""Thermo group additivity estimation: group(O2s-OsCs) + group(O2s-CsH) + group(O2s-(Cds-Cd)H) + group(O2s-OsH) + group(Cs-(Cds-Cds)CsHH) + group(Cs-CsOsOsH) + group(Cds-CdsCsH) + group(Cds-CdsOsH) + radical(Cds_S) + radical(CCOJ) + radical(Cs_P)"""),
)

species(
    label = '[O]O[C]([O])CC=CO(15247)',
    structure = SMILES('[O]O[C]([O])CC=CO'),
    E0 = (69.0883,'kJ/mol'),
    modes = [
        HarmonicOscillator(frequencies=([492.5,1135,1000,360,370,350,2995,3025,975,1000,1300,1375,400,500,1630,1680,3615,1277.5,1000,2750,2850,1437.5,1250,1305,750,350,241.439,241.44,241.443],'cm^-1')),
        HinderedRotor(inertia=(0.0028919,'amu*angstrom^2'), symmetry=1, barrier=(0.119627,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.350135,'amu*angstrom^2'), symmetry=1, barrier=(14.4832,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.35014,'amu*angstrom^2'), symmetry=1, barrier=(14.4832,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.350125,'amu*angstrom^2'), symmetry=1, barrier=(14.4832,'kJ/mol'), semiclassical=False),
    ],
    spinMultiplicity = 4,
    opticalIsomers = 1,
    molecularWeight = (117.08,'amu'),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[0.382292,0.0798171,-9.6622e-05,5.94239e-08,-1.43723e-11,8439.56,31.8599], Tmin=(100,'K'), Tmax=(1012.98,'K')), NASAPolynomial(coeffs=[15.3942,0.0205389,-8.84408e-06,1.65513e-09,-1.15159e-13,5398.21,-40.7574], Tmin=(1012.98,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(69.0883,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(291.007,'J/(mol*K)'), comment="""Thermo group additivity estimation: group(O2s-OsCs) + group(O2s-CsH) + group(O2s-(Cds-Cd)H) + group(O2s-OsH) + group(Cs-(Cds-Cds)CsHH) + group(Cs-CsOsOsH) + group(Cds-CdsCsH) + group(Cds-CdsOsH) + radical(Cs_P) + radical(CCOJ) + radical(ROOJ)"""),
)

species(
    label = '[O]OC(O)[CH][C]=CO(28890)',
    structure = SMILES('[O]OC(O)[CH][C]=CO'),
    E0 = (-7.10477,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
    molecularWeight = (117.08,'amu'),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[-0.525262,0.0935198,-0.000117885,7.08336e-08,-1.62157e-11,-686.323,31.2341], Tmin=(100,'K'), Tmax=(1086.1,'K')), NASAPolynomial(coeffs=[21.9909,0.0105954,-3.35941e-06,5.36554e-10,-3.4696e-14,-5577.32,-79.2533], Tmin=(1086.1,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(-7.10477,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(286.849,'J/(mol*K)'), comment="""Thermo group additivity estimation: group(O2s-OsCs) + group(O2s-CsH) + group(O2s-(Cds-Cd)H) + group(O2s-OsH) + group(Cs-(Cds-Cds)CsHH) + group(Cs-CsOsOsH) + group(Cds-CdsCsH) + group(Cds-CdsOsH) + radical(C=CCJCO) + radical(Cds_S) + radical(ROOJ)"""),
)

species(
    label = '[O]C=CCC([O])O[O](12833)',
    structure = SMILES('[O]C=CCC([O])O[O]'),
    E0 = (5.30452,'kJ/mol'),
    modes = [
        HarmonicOscillator(frequencies=([492.5,1135,1000,2995,3025,975,1000,1300,1375,400,500,1630,1680,1380,1390,370,380,2900,435,2750,2850,1437.5,1250,1305,750,350,331.702,331.703,331.704,331.705],'cm^-1')),
        HinderedRotor(inertia=(0.0218491,'amu*angstrom^2'), symmetry=1, barrier=(1.70595,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.195828,'amu*angstrom^2'), symmetry=1, barrier=(15.2897,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.170128,'amu*angstrom^2'), symmetry=1, barrier=(13.2831,'kJ/mol'), semiclassical=False),
    ],
    spinMultiplicity = 4,
    opticalIsomers = 1,
    molecularWeight = (117.08,'amu'),
    collisionModel = TransportData(shapeIndex=2, epsilon=(4851.45,'J/mol'), sigma=(7.61033,'angstroms'), dipoleMoment=(0,'C*m'), polarizability=(0,'angstroms^3'), rotrelaxcollnum=0, comment="""Epsilon & sigma estimated with Tc=757.79 K, Pc=24.98 bar (from Joback method)"""),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[0.932472,0.0699779,-7.50118e-05,4.29122e-08,-1.00208e-11,746.406,30.6117], Tmin=(100,'K'), Tmax=(1027.07,'K')), NASAPolynomial(coeffs=[11.8885,0.0273084,-1.26935e-05,2.46112e-09,-1.74405e-13,-1504.07,-22.5371], Tmin=(1027.07,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(5.30452,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(295.164,'J/(mol*K)'), comment="""Thermo group additivity estimation: group(O2s-OsCs) + group(O2s-CsH) + group(O2s-(Cds-Cd)H) + group(O2s-OsH) + group(Cs-(Cds-Cds)CsHH) + group(Cs-CsOsOsH) + group(Cds-CdsCsH) + group(Cds-CdsOsH) + radical(CCOJ) + radical(C=COJ) + radical(ROOJ)"""),
)

species(
    label = '[O]C([CH][C]=CO)OO(28891)',
    structure = SMILES('[O]C([CH][C]=CO)OO'),
    E0 = (66.5957,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
    molecularWeight = (117.08,'amu'),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[-0.115109,0.0905574,-0.000111528,6.75671e-08,-1.59837e-11,8157.9,29.8726], Tmin=(100,'K'), Tmax=(1037.14,'K')), NASAPolynomial(coeffs=[18.2336,0.0197922,-9.18308e-06,1.78202e-09,-1.26699e-13,4351.78,-59.3189], Tmin=(1037.14,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(66.5957,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(286.849,'J/(mol*K)'), comment="""Thermo group additivity estimation: group(O2s-OsCs) + group(O2s-CsH) + group(O2s-(Cds-Cd)H) + group(O2s-OsH) + group(Cs-(Cds-Cds)CsHH) + group(Cs-CsOsOsH) + group(Cds-CdsCsH) + group(Cds-CdsOsH) + radical(CCOJ) + radical(C=CCJCO) + radical(Cds_S)"""),
)

species(
    label = '[O]OC(O)C[C]=[C]O(28892)',
    structure = SMILES('[O]OC(O)C[C]=[C]O'),
    E0 = (115.723,'kJ/mol'),
    modes = [
        HarmonicOscillator(frequencies=([1670,1700,300,440,492.5,1135,1000,2750,2850,1437.5,1250,1305,750,350,3580,3650,1210,1345,900,1100,1380,1390,370,380,2900,435,289.804,289.804],'cm^-1')),
        HinderedRotor(inertia=(0.184684,'amu*angstrom^2'), symmetry=1, barrier=(11.0069,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.184684,'amu*angstrom^2'), symmetry=1, barrier=(11.0069,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.184684,'amu*angstrom^2'), symmetry=1, barrier=(11.0069,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.184684,'amu*angstrom^2'), symmetry=1, barrier=(11.0069,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.184684,'amu*angstrom^2'), symmetry=1, barrier=(11.0069,'kJ/mol'), semiclassical=False),
    ],
    spinMultiplicity = 4,
    opticalIsomers = 1,
    molecularWeight = (117.08,'amu'),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[0.0967773,0.0889562,-0.000128405,9.42097e-08,-2.69715e-11,14055.9,34.9371], Tmin=(100,'K'), Tmax=(862.137,'K')), NASAPolynomial(coeffs=[14.9419,0.0200824,-8.57834e-06,1.55381e-09,-1.04279e-13,11496.1,-34.4805], Tmin=(862.137,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(115.723,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(286.849,'J/(mol*K)'), comment="""Thermo group additivity estimation: group(O2s-OsCs) + group(O2s-CsH) + group(O2s-(Cds-Cd)H) + group(O2s-OsH) + group(Cs-(Cds-Cds)CsHH) + group(Cs-CsOsOsH) + group(Cds-CdsCsH) + group(Cds-CdsOsH) + radical(Cds_S) + radical(C=CJO) + radical(ROOJ)"""),
)

species(
    label = '[O]C(C[C]=[C]O)OO(28893)',
    structure = SMILES('[O]C(C[C]=[C]O)OO'),
    E0 = (189.423,'kJ/mol'),
    modes = [
        HarmonicOscillator(frequencies=([1670,1700,300,440,3615,1310,387.5,850,1000,1380,1390,370,380,2900,435,3615,1277.5,1000,2750,2850,1437.5,1250,1305,750,350,200,800,1600],'cm^-1')),
        HinderedRotor(inertia=(0.156089,'amu*angstrom^2'), symmetry=1, barrier=(3.5888,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.156089,'amu*angstrom^2'), symmetry=1, barrier=(3.5888,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.156089,'amu*angstrom^2'), symmetry=1, barrier=(3.5888,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.156089,'amu*angstrom^2'), symmetry=1, barrier=(3.5888,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.156089,'amu*angstrom^2'), symmetry=1, barrier=(3.5888,'kJ/mol'), semiclassical=False),
    ],
    spinMultiplicity = 4,
    opticalIsomers = 1,
    molecularWeight = (117.08,'amu'),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[0.195098,0.0898205,-0.000136216,1.10666e-07,-3.58778e-11,22913.5,34.6849], Tmin=(100,'K'), Tmax=(782.779,'K')), NASAPolynomial(coeffs=[11.797,0.028236,-1.38003e-05,2.6572e-09,-1.84503e-13,21167.6,-17.9965], Tmin=(782.779,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(189.423,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(286.849,'J/(mol*K)'), comment="""Thermo group additivity estimation: group(O2s-OsCs) + group(O2s-CsH) + group(O2s-(Cds-Cd)H) + group(O2s-OsH) + group(Cs-(Cds-Cds)CsHH) + group(Cs-CsOsOsH) + group(Cds-CdsCsH) + group(Cds-CdsOsH) + radical(C=CJO) + radical(CCOJ) + radical(Cds_S)"""),
)

species(
    label = '[O]C=[C]CC(O)O[O](15242)',
    structure = SMILES('[O]C=[C]CC(O)O[O]'),
    E0 = (17.4412,'kJ/mol'),
    modes = [
        HarmonicOscillator(frequencies=([492.5,1135,1000,3615,1277.5,1000,2750,2850,1437.5,1250,1305,750,350,1380,1390,370,380,2900,435,1685,370,3010,987.5,1337.5,450,1655,330.496,333.612,334.548],'cm^-1')),
        HinderedRotor(inertia=(0.16689,'amu*angstrom^2'), symmetry=1, barrier=(13.1239,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.164956,'amu*angstrom^2'), symmetry=1, barrier=(13.1294,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.167,'amu*angstrom^2'), symmetry=1, barrier=(13.1243,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.167225,'amu*angstrom^2'), symmetry=1, barrier=(13.1234,'kJ/mol'), semiclassical=False),
    ],
    spinMultiplicity = 4,
    opticalIsomers = 1,
    molecularWeight = (117.08,'amu'),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[0.375123,0.0793233,-9.49333e-05,5.74251e-08,-1.36319e-11,2228.66,32.2195], Tmin=(100,'K'), Tmax=(1032.6,'K')), NASAPolynomial(coeffs=[15.7486,0.0197702,-8.42332e-06,1.57218e-09,-1.09368e-13,-946.254,-42.4419], Tmin=(1032.6,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(17.4412,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(291.007,'J/(mol*K)'), comment="""Thermo group additivity estimation: group(O2s-OsCs) + group(O2s-CsH) + group(O2s-(Cds-Cd)H) + group(O2s-OsH) + group(Cs-(Cds-Cds)CsHH) + group(Cs-CsOsOsH) + group(Cds-CdsCsH) + group(Cds-CdsOsH) + radical(C=COJ) + radical(ROOJ) + radical(Cds_S)"""),
)

species(
    label = '[O]C=[C]CC([O])OO(15246)',
    structure = SMILES('[O]C=[C]CC([O])OO'),
    E0 = (91.1416,'kJ/mol'),
    modes = [
        HarmonicOscillator(frequencies=([1685,370,3615,1310,387.5,850,1000,3010,987.5,1337.5,450,1655,1380,1390,370,380,2900,435,2750,2850,1437.5,1250,1305,750,350,414.544,415.085,415.67,420.096],'cm^-1')),
        HinderedRotor(inertia=(0.0854338,'amu*angstrom^2'), symmetry=1, barrier=(10.4434,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.0846573,'amu*angstrom^2'), symmetry=1, barrier=(10.4738,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.195087,'amu*angstrom^2'), symmetry=1, barrier=(24.0979,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.302429,'amu*angstrom^2'), symmetry=1, barrier=(37.0275,'kJ/mol'), semiclassical=False),
    ],
    spinMultiplicity = 4,
    opticalIsomers = 1,
    molecularWeight = (117.08,'amu'),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[0.695334,0.0773605,-9.18056e-05,5.80245e-08,-1.49353e-11,11076.9,31.1849], Tmin=(100,'K'), Tmax=(936.748,'K')), NASAPolynomial(coeffs=[12.0406,0.0289152,-1.42312e-05,2.81645e-09,-2.01433e-13,8951.33,-22.8084], Tmin=(936.748,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(91.1416,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(291.007,'J/(mol*K)'), comment="""Thermo group additivity estimation: group(O2s-OsCs) + group(O2s-CsH) + group(O2s-(Cds-Cd)H) + group(O2s-OsH) + group(Cs-(Cds-Cds)CsHH) + group(Cs-CsOsOsH) + group(Cds-CdsCsH) + group(Cds-CdsOsH) + radical(CCOJ) + radical(Cds_S) + radical(C=COJ)"""),
)

species(
    label = '[O][CH]C[C]=CO(28093)',
    structure = SMILES('[O][CH]C[C]=CO'),
    E0 = (258.703,'kJ/mol'),
    modes = [
        HarmonicOscillator(frequencies=([1685,370,2750,2850,1437.5,1250,1305,750,350,3010,987.5,1337.5,450,1655,3615,1277.5,1000,3025,407.5,1350,352.5,305.208,305.211,305.239],'cm^-1')),
        HinderedRotor(inertia=(0.213679,'amu*angstrom^2'), symmetry=1, barrier=(14.1255,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.213686,'amu*angstrom^2'), symmetry=1, barrier=(14.1255,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.21367,'amu*angstrom^2'), symmetry=1, barrier=(14.1253,'kJ/mol'), semiclassical=False),
    ],
    spinMultiplicity = 4,
    opticalIsomers = 1,
    molecularWeight = (85.0813,'amu'),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[1.03666,0.0652876,-8.33377e-05,5.46181e-08,-1.39965e-11,31221.4,24.1632], Tmin=(100,'K'), Tmax=(961.256,'K')), NASAPolynomial(coeffs=[12.9274,0.0158092,-6.1313e-06,1.07434e-09,-7.15047e-14,28935.3,-32.7332], Tmin=(961.256,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(258.703,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(245.277,'J/(mol*K)'), comment="""Thermo group additivity estimation: group(O2s-CsH) + group(O2s-(Cds-Cd)H) + group(Cs-(Cds-Cds)CsHH) + group(Cs-CsOsHH) + group(Cds-CdsCsH) + group(Cds-CdsOsH) + radical(CCsJOH) + radical(Cds_S) + radical(CCOJ)"""),
)

species(
    label = '[CH]=[C]CC([O])O[O](18047)',
    structure = SMILES('[CH]=[C]CC([O])O[O]'),
    E0 = (557.572,'kJ/mol'),
    modes = [
        HarmonicOscillator(frequencies=([1685,370,492.5,1135,1000,1380,1390,370,380,2900,435,3120,650,792.5,1650,2750,2850,1437.5,1250,1305,750,350,180,1748.8],'cm^-1')),
        HinderedRotor(inertia=(0.34039,'amu*angstrom^2'), symmetry=1, barrier=(7.82623,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.340392,'amu*angstrom^2'), symmetry=1, barrier=(7.82628,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.34042,'amu*angstrom^2'), symmetry=1, barrier=(7.82693,'kJ/mol'), semiclassical=False),
    ],
    spinMultiplicity = 5,
    opticalIsomers = 1,
    molecularWeight = (100.073,'amu'),
    collisionModel = TransportData(shapeIndex=2, epsilon=(4216.54,'J/mol'), sigma=(6.85526,'angstroms'), dipoleMoment=(0,'C*m'), polarizability=(0,'angstroms^3'), rotrelaxcollnum=0, comment="""Epsilon & sigma estimated with Tc=658.61 K, Pc=29.7 bar (from Joback method)"""),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[1.16678,0.0728546,-0.000128772,1.23686e-07,-4.52742e-11,67152.3,29.366], Tmin=(100,'K'), Tmax=(847.544,'K')), NASAPolynomial(coeffs=[4.77622,0.0326488,-1.66065e-05,3.20113e-09,-2.20109e-13,67372.7,17.4592], Tmin=(847.544,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(557.572,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(245.277,'J/(mol*K)'), comment="""Thermo group additivity estimation: group(O2s-OsCs) + group(O2s-CsH) + group(O2s-OsH) + group(Cs-(Cds-Cds)CsHH) + group(Cs-CsOsOsH) + group(Cds-CdsCsH) + group(Cds-CdsHH) + radical(CCOJ) + radical(ROOJ) + radical(Cds_S) + radical(Cds_P)"""),
)

species(
    label = '[O]C=[C]CC([O])O[O](15250)',
    structure = SMILES('[O]C=[C]CC([O])O[O]'),
    E0 = (243.146,'kJ/mol'),
    modes = [
        HarmonicOscillator(frequencies=([1685,370,492.5,1135,1000,3010,987.5,1337.5,450,1655,1380,1390,370,380,2900,435,2750,2850,1437.5,1250,1305,750,350,311.181,311.19,311.19,311.197],'cm^-1')),
        HinderedRotor(inertia=(0.225183,'amu*angstrom^2'), symmetry=1, barrier=(15.4734,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.130728,'amu*angstrom^2'), symmetry=1, barrier=(8.98312,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.00174081,'amu*angstrom^2'), symmetry=1, barrier=(0.119627,'kJ/mol'), semiclassical=False),
    ],
    spinMultiplicity = 5,
    opticalIsomers = 1,
    molecularWeight = (116.072,'amu'),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[0.778909,0.0759851,-0.000104714,7.87122e-08,-2.39273e-11,29355,31.6377], Tmin=(100,'K'), Tmax=(801.833,'K')), NASAPolynomial(coeffs=[10.5949,0.0270164,-1.31051e-05,2.54432e-09,-1.78736e-13,27780.9,-13.5505], Tmin=(801.833,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(243.146,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(270.22,'J/(mol*K)'), comment="""Thermo group additivity estimation: group(O2s-OsCs) + group(O2s-CsH) + group(O2s-(Cds-Cd)H) + group(O2s-OsH) + group(Cs-(Cds-Cds)CsHH) + group(Cs-CsOsOsH) + group(Cds-CdsCsH) + group(Cds-CdsOsH) + radical(ROOJ) + radical(C=COJ) + radical(CCOJ) + radical(Cds_S)"""),
)

species(
    label = '[O]OC([O])[CH][C]=CO(28894)',
    structure = SMILES('[O]OC([O])[CH][C]=CO'),
    E0 = (218.6,'kJ/mol'),
    modes = [
        HarmonicOscillator(frequencies=([492.5,1135,1000,3615,1277.5,1000,3025,407.5,1350,352.5,1380,1390,370,380,2900,435,1685,370,3010,987.5,1337.5,450,1655,180,180,180],'cm^-1')),
        HinderedRotor(inertia=(0.937493,'amu*angstrom^2'), symmetry=1, barrier=(21.5548,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.936646,'amu*angstrom^2'), symmetry=1, barrier=(21.5353,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.937775,'amu*angstrom^2'), symmetry=1, barrier=(21.5613,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(1.71614,'amu*angstrom^2'), symmetry=1, barrier=(39.4575,'kJ/mol'), semiclassical=False),
    ],
    spinMultiplicity = 5,
    opticalIsomers = 1,
    molecularWeight = (116.072,'amu'),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[0.120353,0.0872866,-0.000117313,7.82763e-08,-2.03649e-11,26429.6,29.7877], Tmin=(100,'K'), Tmax=(946.416,'K')), NASAPolynomial(coeffs=[16.3649,0.0186256,-8.48447e-06,1.61176e-09,-1.12506e-13,23355,-47.6874], Tmin=(946.416,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(218.6,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(266.063,'J/(mol*K)'), comment="""Thermo group additivity estimation: group(O2s-OsCs) + group(O2s-CsH) + group(O2s-(Cds-Cd)H) + group(O2s-OsH) + group(Cs-(Cds-Cds)CsHH) + group(Cs-CsOsOsH) + group(Cds-CdsCsH) + group(Cds-CdsOsH) + radical(CCOJ) + radical(ROOJ) + radical(Cds_S) + radical(C=CCJCO)"""),
)

species(
    label = '[O]O[C]([O])C[C]=CO(28895)',
    structure = SMILES('[O]O[C]([O])C[C]=CO'),
    E0 = (306.93,'kJ/mol'),
    modes = [
        HarmonicOscillator(frequencies=([492.5,1135,1000,3615,1277.5,1000,360,370,350,2750,2850,1437.5,1250,1305,750,350,1685,370,3010,987.5,1337.5,450,1655,246.508,246.51,246.511],'cm^-1')),
        HinderedRotor(inertia=(0.00277419,'amu*angstrom^2'), symmetry=1, barrier=(0.119627,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.307483,'amu*angstrom^2'), symmetry=1, barrier=(13.2591,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.307482,'amu*angstrom^2'), symmetry=1, barrier=(13.2591,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.307478,'amu*angstrom^2'), symmetry=1, barrier=(13.2591,'kJ/mol'), semiclassical=False),
    ],
    spinMultiplicity = 5,
    opticalIsomers = 1,
    molecularWeight = (116.072,'amu'),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[0.316265,0.084787,-0.000122643,9.02815e-08,-2.60478e-11,37044.4,32.572], Tmin=(100,'K'), Tmax=(853.758,'K')), NASAPolynomial(coeffs=[14.0914,0.0202472,-9.24878e-06,1.73543e-09,-1.19164e-13,34692.3,-31.7071], Tmin=(853.758,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(306.93,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(266.063,'J/(mol*K)'), comment="""Thermo group additivity estimation: group(O2s-OsCs) + group(O2s-CsH) + group(O2s-(Cds-Cd)H) + group(O2s-OsH) + group(Cs-(Cds-Cds)CsHH) + group(Cs-CsOsOsH) + group(Cds-CdsCsH) + group(Cds-CdsOsH) + radical(CCOJ) + radical(Cs_P) + radical(ROOJ) + radical(Cds_S)"""),
)

species(
    label = '[O]OC([O])C[C]=[C]O(28896)',
    structure = SMILES('[O]OC([O])C[C]=[C]O'),
    E0 = (341.428,'kJ/mol'),
    modes = [
        HarmonicOscillator(frequencies=([1670,1700,300,440,492.5,1135,1000,1380,1390,370,380,2900,435,3615,1277.5,1000,2750,2850,1437.5,1250,1305,750,350,328.455,328.457,328.462],'cm^-1')),
        HinderedRotor(inertia=(0.105888,'amu*angstrom^2'), symmetry=1, barrier=(8.10702,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.105898,'amu*angstrom^2'), symmetry=1, barrier=(8.10702,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.105894,'amu*angstrom^2'), symmetry=1, barrier=(8.10697,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.105875,'amu*angstrom^2'), symmetry=1, barrier=(8.10692,'kJ/mol'), semiclassical=False),
    ],
    spinMultiplicity = 5,
    opticalIsomers = 1,
    molecularWeight = (116.072,'amu'),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[0.397409,0.0869976,-0.000143773,1.23802e-07,-4.12661e-11,41186.6,34.7145], Tmin=(100,'K'), Tmax=(855.559,'K')), NASAPolynomial(coeffs=[10.3698,0.0262873,-1.26369e-05,2.37468e-09,-1.60841e-13,39995.7,-8.82804], Tmin=(855.559,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(341.428,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(266.063,'J/(mol*K)'), comment="""Thermo group additivity estimation: group(O2s-OsCs) + group(O2s-CsH) + group(O2s-(Cds-Cd)H) + group(O2s-OsH) + group(Cs-(Cds-Cds)CsHH) + group(Cs-CsOsOsH) + group(Cds-CdsCsH) + group(Cds-CdsOsH) + radical(Cds_S) + radical(ROOJ) + radical(CCOJ) + radical(C=CJO)"""),
)

species(
    label = '[O]OC(O)C=C=CO(28897)',
    structure = SMILES('[O]OC(O)C=C=CO'),
    E0 = (-198.639,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
    molecularWeight = (117.08,'amu'),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[-0.180219,0.0876889,-0.000110314,6.72417e-08,-1.57009e-11,-23736.5,30.2709], Tmin=(100,'K'), Tmax=(1062.49,'K')), NASAPolynomial(coeffs=[19.7699,0.0125823,-4.28086e-06,7.11054e-10,-4.66342e-14,-27975.9,-67.1865], Tmin=(1062.49,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(-198.639,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(291.007,'J/(mol*K)'), comment="""Thermo group additivity estimation: group(O2s-OsCs) + group(O2s-CsH) + group(O2s-(Cds-Cd)H) + group(O2s-OsH) + group(Cs-CsOsOsH) + group(Cds-CdsCsH) + group(Cds-CdsOsH) + group(Cdd-CdsCds) + radical(ROOJ)"""),
)

species(
    label = 'O=C(C[C]=CO)OO(27761)',
    structure = SMILES('O=C(C[C]=CO)OO'),
    E0 = (-331.937,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
    molecularWeight = (117.08,'amu'),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[0.307692,0.0603133,-6.87846e-06,-4.96193e-08,2.66465e-11,-39771,28.7891], Tmin=(100,'K'), Tmax=(969.566,'K')), NASAPolynomial(coeffs=[25.3705,0.00709828,-2.18753e-06,5.38091e-10,-5.10182e-14,-46989.7,-103.514], Tmin=(969.566,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(-331.937,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(286.849,'J/(mol*K)'), comment="""Thermo group additivity estimation: group(O2s-O2s(Cds-O2d)) + group(O2s-(Cds-Cd)H) + group(O2s-OsH) + group(Cs-(Cds-O2d)(Cds-Cds)HH) + group(Cds-CdsCsH) + group(Cds-OdCsOs) + group(Cds-CdsOsH) + radical(Cds_S)"""),
)

species(
    label = '[O]OC(=O)CC=CO(15253)',
    structure = SMILES('[O]OC(=O)CC=CO'),
    E0 = (-375.39,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
    molecularWeight = (117.08,'amu'),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[0.624377,0.0518981,1.39366e-05,-7.04507e-08,3.41882e-11,-45007.1,27.9372], Tmin=(100,'K'), Tmax=(956.861,'K')), NASAPolynomial(coeffs=[24.696,0.0064098,-1.19259e-06,3.1331e-10,-3.49683e-14,-52138,-100.323], Tmin=(956.861,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(-375.39,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(291.007,'J/(mol*K)'), comment="""Thermo group additivity estimation: group(O2s-O2s(Cds-O2d)) + group(O2s-(Cds-Cd)H) + group(O2s-OsH) + group(Cs-(Cds-O2d)(Cds-Cds)HH) + group(Cds-CdsCsH) + group(Cds-OdCsOs) + group(Cds-CdsOsH) + radical(C(=O)OOJ)"""),
)

species(
    label = '[O]C(C=C=CO)OO(28898)',
    structure = SMILES('[O]C(C=C=CO)OO'),
    E0 = (-124.939,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
    molecularWeight = (117.08,'amu'),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[0.203502,0.0850186,-0.000104893,6.50823e-08,-1.59019e-11,-14891.1,29.0056], Tmin=(100,'K'), Tmax=(1001.21,'K')), NASAPolynomial(coeffs=[16.0426,0.0217379,-1.0085e-05,1.95267e-09,-1.38366e-13,-18062.7,-47.4279], Tmin=(1001.21,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(-124.939,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(291.007,'J/(mol*K)'), comment="""Thermo group additivity estimation: group(O2s-OsCs) + group(O2s-CsH) + group(O2s-(Cds-Cd)H) + group(O2s-OsH) + group(Cs-CsOsOsH) + group(Cds-CdsCsH) + group(Cds-CdsOsH) + group(Cdd-CdsCds) + radical(CCOJ)"""),
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
    label = 'OC=[C]CC1OO1(28899)',
    structure = SMILES('OC=[C]CC1OO1'),
    E0 = (46.6799,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
    molecularWeight = (101.081,'amu'),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[0.887249,0.0494094,7.34619e-06,-6.7856e-08,3.68931e-11,5744.36,24.9981], Tmin=(100,'K'), Tmax=(895.387,'K')), NASAPolynomial(coeffs=[24.0431,-0.00175214,5.46653e-06,-1.24218e-09,8.54654e-14,-498.159,-95.8603], Tmin=(895.387,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(46.6799,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(270.22,'J/(mol*K)'), comment="""Thermo group additivity estimation: group(O2s-OsCs) + group(O2s-OsCs) + group(O2s-(Cds-Cd)H) + group(Cs-CsOsOsH) + group(Cs-(Cds-Cds)CsHH) + group(Cds-CdsCsH) + group(Cds-CdsOsH) + ring(dioxirane) + radical(Cds_S)"""),
)

species(
    label = '[O]C1CC(=CO)O1(28900)',
    structure = SMILES('[O]C1CC(=CO)O1'),
    E0 = (-238.556,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
    molecularWeight = (101.081,'amu'),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[1.10963,0.0415519,3.19748e-05,-9.45318e-08,4.67183e-11,-28566.8,19.8656], Tmin=(100,'K'), Tmax=(896.734,'K')), NASAPolynomial(coeffs=[24.4223,-0.00282276,6.48248e-06,-1.44457e-09,9.84393e-14,-35144.8,-103.428], Tmin=(896.734,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(-238.556,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(278.535,'J/(mol*K)'), comment="""Thermo group additivity estimation: group(O2s-Cs(Cds-Cd)) + group(O2s-CsH) + group(O2s-(Cds-Cd)H) + group(Cs-(Cds-Cds)CsHH) + group(Cs-CsOsOsH) + group(Cds-CdsCsOs) + group(Cds-CdsOsH) + ring(2methyleneoxetane) + radical(CCOJ)"""),
)

species(
    label = '[CH2]C(=CO)C([O])O[O](15193)',
    structure = SMILES('[CH2]C(=CO)C([O])O[O]'),
    E0 = (-1.06846,'kJ/mol'),
    modes = [
        HarmonicOscillator(frequencies=([492.5,1135,1000,3615,1277.5,1000,350,440,435,1725,3000,3100,440,815,1455,1000,3010,987.5,1337.5,450,1655,1380,1390,370,380,2900,435,286.812,286.872],'cm^-1')),
        HinderedRotor(inertia=(0.247311,'amu*angstrom^2'), symmetry=1, barrier=(14.4253,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.24731,'amu*angstrom^2'), symmetry=1, barrier=(14.4243,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.247423,'amu*angstrom^2'), symmetry=1, barrier=(14.426,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.555957,'amu*angstrom^2'), symmetry=1, barrier=(32.3788,'kJ/mol'), semiclassical=False),
    ],
    spinMultiplicity = 4,
    opticalIsomers = 1,
    molecularWeight = (117.08,'amu'),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[0.152761,0.0846862,-0.00010561,6.62503e-08,-1.62506e-11,10.0675,30.1205], Tmin=(100,'K'), Tmax=(1001.96,'K')), NASAPolynomial(coeffs=[16.4219,0.0197358,-8.37366e-06,1.55156e-09,-1.07277e-13,-3250.08,-48.4002], Tmin=(1001.96,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(-1.06846,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(291.007,'J/(mol*K)'), comment="""Thermo group additivity estimation: group(O2s-OsCs) + group(O2s-CsH) + group(O2s-(Cds-Cd)H) + group(O2s-OsH) + group(Cs-CsOsOsH) + group(Cs-(Cds-Cds)HHH) + group(Cds-CdsCsCs) + group(Cds-CdsOsH) + radical(CCOJ) + radical(Allyl_P) + radical(ROOJ)"""),
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
    label = '[O]C=C[C]=CO(23190)',
    structure = SMILES('[O]C=C[C]=CO'),
    E0 = (16.794,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
    molecularWeight = (84.0734,'amu'),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[0.93068,0.0480724,2.96789e-06,-6.71246e-08,3.83779e-11,2148.87,19.742], Tmin=(100,'K'), Tmax=(884.221,'K')), NASAPolynomial(coeffs=[26.0788,-0.010819,9.78516e-06,-2.08098e-09,1.44498e-13,-4443.54,-110.619], Tmin=(884.221,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(16.794,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(224.491,'J/(mol*K)'), comment="""Thermo group additivity estimation: group(O2s-(Cds-Cd)H) + group(O2s-(Cds-Cd)H) + group(Cds-Cds(Cds-Cds)H) + group(Cds-Cds(Cds-Cds)H) + group(Cds-CdsOsH) + group(Cds-CdsOsH) + radical(C=CJC=C) + radical(C=COJ)"""),
)

species(
    label = 'OC=[C]CC1OOO1(28901)',
    structure = SMILES('OC=[C]CC1OOO1'),
    E0 = (87.8923,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
    molecularWeight = (117.08,'amu'),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[0.789601,0.0455318,3.51805e-05,-9.58387e-08,4.43984e-11,10709.6,28.441], Tmin=(100,'K'), Tmax=(941.425,'K')), NASAPolynomial(coeffs=[25.5917,0.00354141,1.08272e-06,-1.68057e-10,-1.11148e-15,3230.65,-104.637], Tmin=(941.425,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(87.8923,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(295.164,'J/(mol*K)'), comment="""Thermo group additivity estimation: group(O2s-OsCs) + group(O2s-OsCs) + group(O2s-OsOs) + group(O2s-(Cds-Cd)H) + group(Cs-CsOsOsH) + group(Cs-(Cds-Cds)CsHH) + group(Cds-CdsCsH) + group(Cds-CdsOsH) + ring(Cyclobutane) + radical(Cds_S)"""),
)

species(
    label = '[O]OC1CC(=CO)O1(28874)',
    structure = SMILES('[O]OC1CC(=CO)O1'),
    E0 = (-240.751,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
    molecularWeight = (117.08,'amu'),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[0.338605,0.0575892,8.62066e-06,-8.22727e-08,4.52452e-11,-28802.1,23.2941], Tmin=(100,'K'), Tmax=(890.474,'K')), NASAPolynomial(coeffs=[28.6646,-0.00598436,8.46427e-06,-1.86455e-09,1.29035e-13,-36371,-124.25], Tmin=(890.474,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(-240.751,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(299.321,'J/(mol*K)'), comment="""Thermo group additivity estimation: group(O2s-Cs(Cds-Cd)) + group(O2s-OsCs) + group(O2s-(Cds-Cd)H) + group(O2s-OsH) + group(Cs-(Cds-Cds)CsHH) + group(Cs-CsOsOsH) + group(Cds-CdsCsOs) + group(Cds-CdsOsH) + ring(2methyleneoxetane) + radical(ROOJ)"""),
)

species(
    label = '[O]C1CC(=CO)OO1(28902)',
    structure = SMILES('[O]C1CC(=CO)OO1'),
    E0 = (-203.557,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
    molecularWeight = (117.08,'amu'),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[1.26767,0.0406386,2.82082e-05,-7.38124e-08,3.29949e-11,-24366.2,27.2217], Tmin=(100,'K'), Tmax=(959.393,'K')), NASAPolynomial(coeffs=[19.4249,0.0130542,-3.89713e-06,7.75322e-10,-6.39958e-14,-30064.7,-71.1652], Tmin=(959.393,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(-203.557,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(303.478,'J/(mol*K)'), comment="""Thermo group additivity estimation: group(O2s-OsCs) + group(O2s-O2s(Cds-Cd)) + group(O2s-CsH) + group(O2s-(Cds-Cd)H) + group(Cs-(Cds-Cds)CsHH) + group(Cs-CsOsOsH) + group(Cds-CdsCsOs) + group(Cds-CdsOsH) + ring(Cyclopentane) + radical(CCOJ)"""),
)

species(
    label = '[O]C([O])C[C]=CO(28903)',
    structure = SMILES('[O]C([O])C[C]=CO'),
    E0 = (103.879,'kJ/mol'),
    modes = [
        HarmonicOscillator(frequencies=([1685,370,2750,2850,1437.5,1250,1305,750,350,3010,987.5,1337.5,450,1655,3615,1277.5,1000,1380,1390,370,380,2900,435,446.058,446.063,446.09,446.156],'cm^-1')),
        HinderedRotor(inertia=(0.0681442,'amu*angstrom^2'), symmetry=1, barrier=(9.62277,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.0681446,'amu*angstrom^2'), symmetry=1, barrier=(9.62249,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.068145,'amu*angstrom^2'), symmetry=1, barrier=(9.6229,'kJ/mol'), semiclassical=False),
    ],
    spinMultiplicity = 4,
    opticalIsomers = 1,
    molecularWeight = (101.081,'amu'),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[1.15753,0.0666507,-8.35841e-05,5.82402e-08,-1.65898e-11,12592.5,27.517], Tmin=(100,'K'), Tmax=(850.371,'K')), NASAPolynomial(coeffs=[9.74613,0.0262524,-1.23257e-05,2.37695e-09,-1.67027e-13,11131.7,-12.5262], Tmin=(850.371,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(103.879,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(270.22,'J/(mol*K)'), comment="""Thermo group additivity estimation: group(O2s-CsH) + group(O2s-CsH) + group(O2s-(Cds-Cd)H) + group(Cs-(Cds-Cds)CsHH) + group(Cs-CsOsOsH) + group(Cds-CdsCsH) + group(Cds-CdsOsH) + radical(CCOJ) + radical(CCOJ) + radical(Cds_S)"""),
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
    label = '[C]=CO(27807)',
    structure = SMILES('[C]=CO'),
    E0 = (391.459,'kJ/mol'),
    modes = [
        HarmonicOscillator(frequencies=([3615,1277.5,1000,3010,987.5,1337.5,450,1655],'cm^-1')),
        HinderedRotor(inertia=(1.23325,'amu*angstrom^2'), symmetry=1, barrier=(28.3547,'kJ/mol'), semiclassical=False),
    ],
    spinMultiplicity = 3,
    opticalIsomers = 1,
    molecularWeight = (42.0367,'amu'),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[2.88904,0.0147275,1.32236e-05,-4.12494e-08,2.11475e-11,47130.8,9.58665], Tmin=(100,'K'), Tmax=(884.362,'K')), NASAPolynomial(coeffs=[13.839,-0.00784764,5.80008e-06,-1.19218e-09,8.19757e-14,44140.2,-47.8537], Tmin=(884.362,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(391.459,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(103.931,'J/(mol*K)'), comment="""Thermo group additivity estimation: group(O2s-(Cds-Cd)H) + group(Cds-CdsOsH) + group(Cds-CdsHH) + radical(CdCdJ2_triplet)"""),
)

species(
    label = '[O]O[CH]C[C]=CO(28904)',
    structure = SMILES('[O]O[CH]C[C]=CO'),
    E0 = (264.792,'kJ/mol'),
    modes = [
        HarmonicOscillator(frequencies=([492.5,1135,1000,3615,1277.5,1000,3025,407.5,1350,352.5,2750,2850,1437.5,1250,1305,750,350,1685,370,3010,987.5,1337.5,450,1655,180,180],'cm^-1')),
        HinderedRotor(inertia=(0.727342,'amu*angstrom^2'), symmetry=1, barrier=(16.723,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.727318,'amu*angstrom^2'), symmetry=1, barrier=(16.7225,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.72731,'amu*angstrom^2'), symmetry=1, barrier=(16.7223,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.727293,'amu*angstrom^2'), symmetry=1, barrier=(16.7219,'kJ/mol'), semiclassical=False),
    ],
    spinMultiplicity = 4,
    opticalIsomers = 1,
    molecularWeight = (101.081,'amu'),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[0.416444,0.076827,-9.78785e-05,6.22582e-08,-1.52959e-11,31977.9,27.6464], Tmin=(100,'K'), Tmax=(1007.5,'K')), NASAPolynomial(coeffs=[16.1215,0.0144744,-5.04593e-06,8.3064e-10,-5.33466e-14,28813.4,-48.2388], Tmin=(1007.5,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(264.792,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(266.063,'J/(mol*K)'), comment="""Thermo group additivity estimation: group(O2s-OsCs) + group(O2s-(Cds-Cd)H) + group(O2s-OsH) + group(Cs-(Cds-Cds)CsHH) + group(Cs-CsOsHH) + group(Cds-CdsCsH) + group(Cds-CdsOsH) + radical(ROOJ) + radical(CCsJOOH) + radical(Cds_S)"""),
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
    E0 = (101.684,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS2',
    E0 = (247.364,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS3',
    E0 = (120.762,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS4',
    E0 = (324.187,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS5',
    E0 = (219.434,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS6',
    E0 = (101.684,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS7',
    E0 = (267.132,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS8',
    E0 = (101.684,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS9',
    E0 = (259.363,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS10',
    E0 = (259.72,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS11',
    E0 = (335.741,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS12',
    E0 = (290.097,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS13',
    E0 = (243.661,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS14',
    E0 = (176.964,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS15',
    E0 = (252.029,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS16',
    E0 = (231.855,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS17',
    E0 = (148.763,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS18',
    E0 = (366.665,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS19',
    E0 = (174.64,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS20',
    E0 = (178.072,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS21',
    E0 = (250.076,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS22',
    E0 = (586.072,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS23',
    E0 = (454.951,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS24',
    E0 = (415.56,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS25',
    E0 = (430.405,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS26',
    E0 = (518.735,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS27',
    E0 = (553.233,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS28',
    E0 = (179.931,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS29',
    E0 = (165.084,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS30',
    E0 = (190.652,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS31',
    E0 = (119.466,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS32',
    E0 = (101.684,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS33',
    E0 = (291.214,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS34',
    E0 = (184.582,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS35',
    E0 = (271.721,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS36',
    E0 = (241.053,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS37',
    E0 = (109.968,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS38',
    E0 = (109.968,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS39',
    E0 = (108.796,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS40',
    E0 = (352.08,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS41',
    E0 = (632.15,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS42',
    E0 = (671.673,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

reaction(
    label = 'reaction1',
    reactants = ['[O]OC([O])C[C]=CO(15241)'],
    products = ['[O]OC=O(5472)', 'C=C=CO(12571)'],
    transitionState = 'TS1',
    kinetics = Arrhenius(A=(5e+12,'s^-1'), n=0, Ea=(0,'kJ/mol'), T0=(1,'K'), Tmin=(300,'K'), Tmax=(1500,'K'), comment="""Exact match found for rate rule [RJJ]
Euclidian distance = 0
family: 1,4_Linear_birad_scission"""),
)

reaction(
    label = 'reaction2',
    reactants = ['H(8)', '[O]OC([O])C=C=CO(28885)'],
    products = ['[O]OC([O])C[C]=CO(15241)'],
    transitionState = 'TS2',
    kinetics = Arrhenius(A=(1.23e+08,'cm^3/(mol*s)'), n=1.64, Ea=(8.49352,'kJ/mol'), T0=(1,'K'), Tmin=(300,'K'), Tmax=(1500,'K'), comment="""From training reaction 2566 used for Cds-CsH_Ca;HJ
Exact match found for rate rule [Cds-CsH_Ca;HJ]
Euclidian distance = 0
family: R_Addition_MultipleBond"""),
)

reaction(
    label = 'reaction3',
    reactants = ['H(8)', '[O]OC(=O)C[C]=CO(28886)'],
    products = ['[O]OC([O])C[C]=CO(15241)'],
    transitionState = 'TS3',
    kinetics = Arrhenius(A=(1.83701,'m^3/(mol*s)'), n=1.71338, Ea=(46.5052,'kJ/mol'), T0=(1,'K'), comment="""Estimated using average of templates [CO_O;HJ] + [CO-NdNd_O;YJ] for rate rule [CO-NdNd_O;HJ]
Euclidian distance = 1.0
family: R_Addition_MultipleBond"""),
)

reaction(
    label = 'reaction4',
    reactants = ['H(8)', '[O]OC([O])CC#CO(28887)'],
    products = ['[O]OC([O])C[C]=CO(15241)'],
    transitionState = 'TS4',
    kinetics = Arrhenius(A=(4278.27,'m^3/(mol*s)'), n=1.383, Ea=(15.1097,'kJ/mol'), T0=(1,'K'), comment="""Estimated using template [Ct_Ct;HJ] for rate rule [Ct-O_Ct;HJ]
Euclidian distance = 1.0
family: R_Addition_MultipleBond"""),
)

reaction(
    label = 'reaction5',
    reactants = ['[O][CH]O[O](8201)', 'C=C=CO(12571)'],
    products = ['[O]OC([O])C[C]=CO(15241)'],
    transitionState = 'TS5',
    kinetics = Arrhenius(A=(0.00429749,'m^3/(mol*s)'), n=2.45395, Ea=(16.6104,'kJ/mol'), T0=(1,'K'), comment="""Estimated using an average for rate rule [Cds-HH_Ca;CJ]
Euclidian distance = 0
family: R_Addition_MultipleBond"""),
)

reaction(
    label = 'reaction6',
    reactants = ['[O]OC=O(5472)', '[CH2][C]=CO(18753)'],
    products = ['[O]OC([O])C[C]=CO(15241)'],
    transitionState = 'TS6',
    kinetics = Arrhenius(A=(6.87291,'m^3/(mol*s)'), n=1.39198, Ea=(110.507,'kJ/mol'), T0=(1,'K'), comment="""Estimated using average of templates [CO_O;CJ] + [CO-NdH_O;YJ] for rate rule [CO-NdH_O;CJ]
Euclidian distance = 1.0
family: R_Addition_MultipleBond
Ea raised from 108.5 to 110.5 kJ/mol to match endothermicity of reaction."""),
)

reaction(
    label = 'reaction7',
    reactants = ['OH(D)(132)', 'C#CCC([O])O[O](18040)'],
    products = ['[O]OC([O])C[C]=CO(15241)'],
    transitionState = 'TS7',
    kinetics = Arrhenius(A=(6.508e+07,'cm^3/(mol*s)'), n=1.628, Ea=(0,'kJ/mol'), T0=(1,'K'), comment="""From training reaction 211 used for Ct-H_Ct-Cs;OJ_pri
Exact match found for rate rule [Ct-H_Ct-Cs;OJ_pri]
Euclidian distance = 0
family: R_Addition_MultipleBond
Ea raised from -1.9 to 0 kJ/mol."""),
)

reaction(
    label = 'reaction8',
    reactants = ['O2(2)', 'O=CC[C]=CO(27745)'],
    products = ['[O]OC([O])C[C]=CO(15241)'],
    transitionState = 'TS8',
    kinetics = Arrhenius(A=(0.0849,'cm^3/(mol*s)'), n=3.486, Ea=(174.4,'kJ/mol'), T0=(1,'K'), Tmin=(300,'K'), Tmax=(1600,'K'), comment="""Estimated using template [CO-CsH_O;OJ] for rate rule [CO-CsH_O;O2b]
Euclidian distance = 1.0
Multiplied by reaction path degeneracy 2.0
family: R_Addition_MultipleBond
Ea raised from 172.4 to 174.4 kJ/mol to match endothermicity of reaction."""),
)

reaction(
    label = 'reaction9',
    reactants = ['[O]OC([O])C[C]=CO(15241)'],
    products = ['[O]OC([O])[CH]C=CO(15244)'],
    transitionState = 'TS9',
    kinetics = Arrhenius(A=(2.66329e+10,'s^-1'), n=0.993, Ea=(157.679,'kJ/mol'), T0=(1,'K'), comment="""Estimated using an average for rate rule [R2H_S;Cd_rad_out_Cd;Cs_H_out_H/NonDeC]
Euclidian distance = 0
Multiplied by reaction path degeneracy 2.0
family: intra_H_migration"""),
)

reaction(
    label = 'reaction10',
    reactants = ['[O]OC([O])C[C]=CO(15241)'],
    products = ['[O]O[C](O)C[C]=CO(28888)'],
    transitionState = 'TS10',
    kinetics = Arrhenius(A=(1.43381e+07,'s^-1'), n=1.70481, Ea=(158.036,'kJ/mol'), T0=(1,'K'), comment="""Estimated using template [R2H_S;Y_rad_out;Cs_H_out_NonDe] for rate rule [R2H_S;O_rad_out;Cs_H_out_NDMustO]
Euclidian distance = 1.41421356237
family: intra_H_migration"""),
)

reaction(
    label = 'reaction11',
    reactants = ['[O]OC([O])C[C]=CO(15241)'],
    products = ['[O]OC([O])CC=[C]O(15239)'],
    transitionState = 'TS11',
    kinetics = Arrhenius(A=(1.231e+11,'s^-1'), n=0.765, Ea=(234.057,'kJ/mol'), T0=(1,'K'), comment="""From training reaction 82 used for R2H_D;Cd_rad_out_Cs;Cd_H_out_singleNd
Exact match found for rate rule [R2H_D;Cd_rad_out_Cs;Cd_H_out_singleNd]
Euclidian distance = 0
family: intra_H_migration"""),
)

reaction(
    label = 'reaction12',
    reactants = ['[O][C](C[C]=CO)OO(28889)'],
    products = ['[O]OC([O])C[C]=CO(15241)'],
    transitionState = 'TS12',
    kinetics = Arrhenius(A=(40813.3,'s^-1'), n=2.17068, Ea=(135.172,'kJ/mol'), T0=(1,'K'), comment="""Estimated using average of templates [R3H_SS;Y_rad_out;O_H_out] + [R3H_SS_O;Y_rad_out;XH_out] for rate rule [R3H_SS_O;Y_rad_out;O_H_out]
Euclidian distance = 1.0
family: intra_H_migration"""),
)

reaction(
    label = 'reaction13',
    reactants = ['[O]OC([O])C[C]=CO(15241)'],
    products = ['[O]O[C]([O])CC=CO(15247)'],
    transitionState = 'TS13',
    kinetics = Arrhenius(A=(2.4115e+09,'s^-1'), n=1.00333, Ea=(141.977,'kJ/mol'), T0=(1,'K'), comment="""Estimated using an average for rate rule [R3H_SS_Cs;Cd_rad_out_Cd;XH_out]
Euclidian distance = 0
family: intra_H_migration"""),
)

reaction(
    label = 'reaction14',
    reactants = ['[O]OC([O])C[C]=CO(15241)'],
    products = ['[O]OC(O)[CH][C]=CO(28890)'],
    transitionState = 'TS14',
    kinetics = Arrhenius(A=(223829,'s^-1'), n=2.27675, Ea=(75.2806,'kJ/mol'), T0=(1,'K'), comment="""Estimated using an average for rate rule [R3H_SS_Cs;O_rad_out;XH_out]
Euclidian distance = 0
Multiplied by reaction path degeneracy 2.0
family: intra_H_migration"""),
)

reaction(
    label = 'reaction12',
    reactants = ['[O]OC([O])C[C]=CO(15241)'],
    products = ['[O]C=CCC([O])O[O](12833)'],
    transitionState = 'TS15',
    kinetics = Arrhenius(A=(4.96975e+09,'s^-1'), n=0.933333, Ea=(150.345,'kJ/mol'), T0=(1,'K'), comment="""Estimated using template [R3H_DS;Cd_rad_out_Cs;XH_out] for rate rule [R3H_DS;Cd_rad_out_Cs;O_H_out]
Euclidian distance = 1.0
family: intra_H_migration"""),
)

reaction(
    label = 'reaction16',
    reactants = ['[O]OC([O])C[C]=CO(15241)'],
    products = ['[O]C([CH][C]=CO)OO(28891)'],
    transitionState = 'TS16',
    kinetics = Arrhenius(A=(2.34273e+07,'s^-1'), n=1.54267, Ea=(130.172,'kJ/mol'), T0=(1,'K'), comment="""Estimated using an average for rate rule [R4H_SSS_OCs;O_rad_out;XH_out]
Euclidian distance = 0
Multiplied by reaction path degeneracy 2.0
family: intra_H_migration"""),
)

reaction(
    label = 'reaction17',
    reactants = ['[O]OC(O)C[C]=[C]O(28892)'],
    products = ['[O]OC([O])C[C]=CO(15241)'],
    transitionState = 'TS17',
    kinetics = Arrhenius(A=(136000,'s^-1'), n=1.9199, Ea=(33.0402,'kJ/mol'), T0=(1,'K'), comment="""Estimated using template [R5Hall;Cd_rad_out_single;XH_out] for rate rule [R5HJ_1;Cd_rad_out_singleNd;O_H_out]
Euclidian distance = 1.73205080757
family: intra_H_migration"""),
)

reaction(
    label = 'reaction18',
    reactants = ['[O]C(C[C]=[C]O)OO(28893)'],
    products = ['[O]OC([O])C[C]=CO(15241)'],
    transitionState = 'TS18',
    kinetics = Arrhenius(A=(8.44313e+09,'s^-1'), n=0.985167, Ea=(177.241,'kJ/mol'), T0=(1,'K'), comment="""Estimated using template [RnH;Cd_rad_out_singleNd;XH_out] for rate rule [R6HJ_1;Cd_rad_out_singleNd;O_H_out]
Euclidian distance = 2.2360679775
family: intra_H_migration"""),
)

reaction(
    label = 'reaction19',
    reactants = ['[O]OC([O])C[C]=CO(15241)'],
    products = ['[O]C=[C]CC(O)O[O](15242)'],
    transitionState = 'TS19',
    kinetics = Arrhenius(A=(1.74437e+06,'s^-1'), n=0.972854, Ea=(72.9565,'kJ/mol'), T0=(1,'K'), comment="""Estimated using template [R6Hall;O_rad_out;XH_out] for rate rule [R6HJ_3;O_rad_out;O_H_out]
Euclidian distance = 1.41421356237
family: intra_H_migration"""),
)

reaction(
    label = 'reaction20',
    reactants = ['[O]OC([O])C[C]=CO(15241)'],
    products = ['[O]C=[C]CC([O])OO(15246)'],
    transitionState = 'TS20',
    kinetics = Arrhenius(A=(970995,'s^-1'), n=0.905106, Ea=(76.3888,'kJ/mol'), T0=(1,'K'), comment="""Estimated using template [R7Hall;O_rad_out;XH_out] for rate rule [R7HJ_4;O_rad_out;O_H_out]
Euclidian distance = 1.41421356237
family: intra_H_migration"""),
)

reaction(
    label = 'reaction21',
    reactants = ['O2(2)', '[O][CH]C[C]=CO(28093)'],
    products = ['[O]OC([O])C[C]=CO(15241)'],
    transitionState = 'TS21',
    kinetics = Arrhenius(A=(2.18266e+06,'m^3/(mol*s)'), n=0.193158, Ea=(0,'kJ/mol'), T0=(1,'K'), comment="""Estimated using an average for rate rule [Y_rad;O2_birad]
Euclidian distance = 0
Multiplied by reaction path degeneracy 2.0
family: R_Recombination
Ea raised from -25.1 to 0 kJ/mol."""),
)

reaction(
    label = 'reaction22',
    reactants = ['OH(D)(132)', '[CH]=[C]CC([O])O[O](18047)'],
    products = ['[O]OC([O])C[C]=CO(15241)'],
    transitionState = 'TS22',
    kinetics = Arrhenius(A=(3.05166e+07,'m^3/(mol*s)'), n=0.045, Ea=(0.1046,'kJ/mol'), T0=(1,'K'), comment="""Estimated using an average for rate rule [O_pri_rad;Y_rad]
Euclidian distance = 0
family: R_Recombination"""),
)

reaction(
    label = 'reaction23',
    reactants = ['H(8)', '[O]C=[C]CC([O])O[O](15250)'],
    products = ['[O]OC([O])C[C]=CO(15241)'],
    transitionState = 'TS23',
    kinetics = Arrhenius(A=(7e+11,'cm^3/(mol*s)'), n=0.493, Ea=(0,'kJ/mol'), T0=(1,'K'), Tmin=(200,'K'), Tmax=(2500,'K'), comment="""From training reaction 25 used for H_rad;O_rad/OneDe
Exact match found for rate rule [H_rad;O_rad/OneDe]
Euclidian distance = 0
family: R_Recombination
Ea raised from -1.2 to 0 kJ/mol."""),
)

reaction(
    label = 'reaction24',
    reactants = ['[O][CH]O[O](8201)', '[CH2][C]=CO(18753)'],
    products = ['[O]OC([O])C[C]=CO(15241)'],
    transitionState = 'TS24',
    kinetics = Arrhenius(A=(1.9789e+07,'m^3/(mol*s)'), n=-0.126319, Ea=(0,'kJ/mol'), T0=(1,'K'), comment="""Estimated using an average for rate rule [Y_rad;Y_rad]
Euclidian distance = 0
family: R_Recombination
Ea raised from -15.6 to -15.6 kJ/mol.
Ea raised from -15.6 to 0 kJ/mol."""),
)

reaction(
    label = 'reaction25',
    reactants = ['H(8)', '[O]OC([O])[CH][C]=CO(28894)'],
    products = ['[O]OC([O])C[C]=CO(15241)'],
    transitionState = 'TS25',
    kinetics = Arrhenius(A=(4.34078e+06,'m^3/(mol*s)'), n=0.278577, Ea=(0,'kJ/mol'), T0=(1,'K'), comment="""Estimated using an average for rate rule [Y_rad;H_rad]
Euclidian distance = 0
family: R_Recombination
Ea raised from -1.4 to 0 kJ/mol."""),
)

reaction(
    label = 'reaction26',
    reactants = ['H(8)', '[O]O[C]([O])C[C]=CO(28895)'],
    products = ['[O]OC([O])C[C]=CO(15241)'],
    transitionState = 'TS26',
    kinetics = Arrhenius(A=(4.34078e+06,'m^3/(mol*s)'), n=0.278577, Ea=(0,'kJ/mol'), T0=(1,'K'), comment="""Estimated using an average for rate rule [Y_rad;H_rad]
Euclidian distance = 0
family: R_Recombination
Ea raised from -1.4 to 0 kJ/mol."""),
)

reaction(
    label = 'reaction27',
    reactants = ['H(8)', '[O]OC([O])C[C]=[C]O(28896)'],
    products = ['[O]OC([O])C[C]=CO(15241)'],
    transitionState = 'TS27',
    kinetics = Arrhenius(A=(4.34078e+06,'m^3/(mol*s)'), n=0.278577, Ea=(0,'kJ/mol'), T0=(1,'K'), comment="""Estimated using an average for rate rule [Y_rad;H_rad]
Euclidian distance = 0
family: R_Recombination
Ea raised from -1.4 to 0 kJ/mol."""),
)

reaction(
    label = 'reaction28',
    reactants = ['[O]OC([O])C[C]=CO(15241)'],
    products = ['[O]OC(O)C=C=CO(28897)'],
    transitionState = 'TS28',
    kinetics = Arrhenius(A=(4.00798e+09,'s^-1'), n=0.37, Ea=(78.2471,'kJ/mol'), T0=(1,'K'), comment="""Estimated using average of templates [R3;Y_rad;XH_Rrad_De] + [R3radExo;Y_rad;XH_Rrad] for rate rule [R3radExo;Y_rad;XH_Rrad_De]
Euclidian distance = 1.0
Multiplied by reaction path degeneracy 2.0
family: Intra_Disproportionation"""),
)

reaction(
    label = 'reaction29',
    reactants = ['[O]OC([O])C[C]=CO(15241)'],
    products = ['O=C(C[C]=CO)OO(27761)'],
    transitionState = 'TS29',
    kinetics = Arrhenius(A=(7.437e+08,'s^-1'), n=1.045, Ea=(63.4002,'kJ/mol'), T0=(1,'K'), comment="""Estimated using an average for rate rule [R3radExo;Y_rad;XH_Rrad]
Euclidian distance = 0
family: Intra_Disproportionation"""),
)

reaction(
    label = 'reaction30',
    reactants = ['[O]OC([O])C[C]=CO(15241)'],
    products = ['[O]OC(=O)CC=CO(15253)'],
    transitionState = 'TS30',
    kinetics = Arrhenius(A=(2.6374e+09,'s^-1'), n=0.37, Ea=(88.9686,'kJ/mol'), T0=(1,'K'), comment="""Estimated using average of templates [R3;Y_rad_De;XH_Rrad] + [R3radExo;Y_rad;XH_Rrad] for rate rule [R3radExo;Y_rad_De;XH_Rrad]
Euclidian distance = 1.0
family: Intra_Disproportionation"""),
)

reaction(
    label = 'reaction31',
    reactants = ['[O]OC([O])C[C]=CO(15241)'],
    products = ['[O]C(C=C=CO)OO(28898)'],
    transitionState = 'TS31',
    kinetics = Arrhenius(A=(5.18e+08,'s^-1'), n=0.311, Ea=(17.782,'kJ/mol'), T0=(1,'K'), Tmin=(300,'K'), Tmax=(1500,'K'), comment="""Estimated using template [R4;Y_rad;XH_Rrad_De] for rate rule [R4radExo;Y_rad;XH_Rrad_De]
Euclidian distance = 1.0
Multiplied by reaction path degeneracy 2.0
family: Intra_Disproportionation"""),
)

reaction(
    label = 'reaction32',
    reactants = ['[O]OC([O])C[C]=CO(15241)'],
    products = ['O2(S)(5486)', 'O=CC[C]=CO(27745)'],
    transitionState = 'TS32',
    kinetics = Arrhenius(A=(5e+12,'s^-1'), n=0, Ea=(0,'kJ/mol'), T0=(1,'K'), Tmin=(300,'K'), Tmax=(1500,'K'), comment="""Exact match found for rate rule [RJJ]
Euclidian distance = 0
family: 1,4_Linear_birad_scission"""),
)

reaction(
    label = 'reaction33',
    reactants = ['[O]OC([O])C[C]=CO(15241)'],
    products = ['O(T)(63)', 'OC=[C]CC1OO1(28899)'],
    transitionState = 'TS33',
    kinetics = Arrhenius(A=(2.30161e+09,'s^-1'), n=1.08, Ea=(189.53,'kJ/mol'), T0=(1,'K'), comment="""Estimated using an average for rate rule [R2OO_S;Y_rad_intra;OOJ]
Euclidian distance = 0
family: Cyclic_Ether_Formation"""),
)

reaction(
    label = 'reaction34',
    reactants = ['[O]OC([O])C[C]=CO(15241)'],
    products = ['O(T)(63)', '[O]C1CC(=CO)O1(28900)'],
    transitionState = 'TS34',
    kinetics = Arrhenius(A=(3.11355e+11,'s^-1'), n=0, Ea=(82.8981,'kJ/mol'), T0=(1,'K'), comment="""Estimated using template [R3OO_SS;Y_rad_intra;OO] for rate rule [R3OO_SS;Cd_rad_out;OOJ]
Euclidian distance = 1.41421356237
family: Cyclic_Ether_Formation"""),
)

reaction(
    label = 'reaction35',
    reactants = ['[O]OC([O])C[C]=CO(15241)'],
    products = ['[CH2]C(=CO)C([O])O[O](15193)'],
    transitionState = 'TS35',
    kinetics = Arrhenius(A=(1.74842e+09,'s^-1'), n=1.084, Ea=(170.038,'kJ/mol'), T0=(1,'K'), comment="""Estimated using average of templates [cCsCJ;CdsJ;C] + [cCs(-HH)CJ;CJ;C] for rate rule [cCs(-HH)CJ;CdsJ;C]
Euclidian distance = 1.0
family: 1,2_shiftC"""),
)

reaction(
    label = 'reaction36',
    reactants = ['[O]OC([O])C[C]=CO(15241)'],
    products = ['[O]O(16)', '[O]C=C[C]=CO(23190)'],
    transitionState = 'TS36',
    kinetics = Arrhenius(A=(1.91635e+11,'s^-1'), n=0.573333, Ea=(139.369,'kJ/mol'), T0=(1,'K'), comment="""Estimated using an average for rate rule [R2OO]
Euclidian distance = 0
Multiplied by reaction path degeneracy 2.0
family: HO2_Elimination_from_PeroxyRadical"""),
)

reaction(
    label = 'reaction37',
    reactants = ['[O]OC([O])C[C]=CO(15241)'],
    products = ['OC=[C]CC1OOO1(28901)'],
    transitionState = 'TS37',
    kinetics = Arrhenius(A=(1.62e+12,'s^-1'), n=-0.305, Ea=(8.28432,'kJ/mol'), T0=(1,'K'), comment="""Estimated using template [R4_SSS;Y_rad_out;Ypri_rad_out] for rate rule [R4_SSS;O_rad;Opri_rad]
Euclidian distance = 1.41421356237
family: Birad_recombination"""),
)

reaction(
    label = 'reaction38',
    reactants = ['[O]OC([O])C[C]=CO(15241)'],
    products = ['[O]OC1CC(=CO)O1(28874)'],
    transitionState = 'TS38',
    kinetics = Arrhenius(A=(1.62e+12,'s^-1'), n=-0.305, Ea=(8.28432,'kJ/mol'), T0=(1,'K'), comment="""Estimated using template [R4_SSS;Y_rad_out;Ypri_rad_out] for rate rule [R4_SSS;Y_rad_out;Opri_rad]
Euclidian distance = 1.0
family: Birad_recombination"""),
)

reaction(
    label = 'reaction39',
    reactants = ['[O]OC([O])C[C]=CO(15241)'],
    products = ['[O]C1CC(=CO)OO1(28902)'],
    transitionState = 'TS39',
    kinetics = Arrhenius(A=(7.76e+09,'s^-1'), n=0.311, Ea=(7.1128,'kJ/mol'), T0=(1,'K'), comment="""Estimated using template [R5_SSSS;Y_rad_out;Ypri_rad_out] for rate rule [R5_SSSS;O_rad;Ypri_rad_out]
Euclidian distance = 1.0
family: Birad_recombination"""),
)

reaction(
    label = 'reaction40',
    reactants = ['O(T)(63)', '[O]C([O])C[C]=CO(28903)'],
    products = ['[O]OC([O])C[C]=CO(15241)'],
    transitionState = 'TS40',
    kinetics = Arrhenius(A=(87.1677,'m^3/(mol*s)'), n=1.88017, Ea=(5.1666,'kJ/mol'), T0=(1,'K'), Tmin=(303.03,'K'), Tmax=(2000,'K'), comment="""From training reaction 3 used for O_rad/NonDe;O_birad
Exact match found for rate rule [O_rad/NonDe;O_birad]
Euclidian distance = 0
Multiplied by reaction path degeneracy 2.0
family: Birad_R_Recombination"""),
)

reaction(
    label = 'reaction41',
    reactants = ['[CH2]C([O])O[O](607)', '[C]=CO(27807)'],
    products = ['[O]OC([O])C[C]=CO(15241)'],
    transitionState = 'TS41',
    kinetics = Arrhenius(A=(1.14854e+06,'m^3/(mol*s)'), n=0.575199, Ea=(34.3157,'kJ/mol'), T0=(1,'K'), comment="""Estimated using template [Y_rad;Birad] for rate rule [C_rad/H2/Cs;Birad]
Euclidian distance = 3.0
family: Birad_R_Recombination"""),
)

reaction(
    label = 'reaction42',
    reactants = ['O(T)(63)', '[O]O[CH]C[C]=CO(28904)'],
    products = ['[O]OC([O])C[C]=CO(15241)'],
    transitionState = 'TS42',
    kinetics = Arrhenius(A=(93609.6,'m^3/(mol*s)'), n=1.13083, Ea=(163.847,'kJ/mol'), T0=(1,'K'), Tmin=(303.03,'K'), Tmax=(2000,'K'), comment="""Estimated using template [Y_rad;O_birad] for rate rule [C_rad/H/CsO;O_birad]
Euclidian distance = 4.0
family: Birad_R_Recombination"""),
)

network(
    label = '5160',
    isomers = [
        '[O]OC([O])C[C]=CO(15241)',
    ],
    reactants = [
        ('[O]OC=O(5472)', 'C=C=CO(12571)'),
    ],
    bathGas = {
        'N2': 0.25,
        'Ne': 0.25,
        'He': 0.25,
        'Ar': 0.25,
    },
)

pressureDependence(
    label = '5160',
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

