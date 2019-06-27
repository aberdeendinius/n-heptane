species(
    label = 'C=[C]C(O)C([O])O[O](15278)',
    structure = SMILES('C=[C]C(O)C([O])O[O]'),
    E0 = (137.67,'kJ/mol'),
    modes = [
        HarmonicOscillator(frequencies=([1685,370,492.5,1135,1000,2950,3100,1380,975,1025,1650,3615,1277.5,1000,1380,1383.33,1386.67,1390,370,373.333,376.667,380,2800,3000,430,440,351.961,351.961,1965.7],'cm^-1')),
        HinderedRotor(inertia=(0.10238,'amu*angstrom^2'), symmetry=1, barrier=(8.99978,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.10238,'amu*angstrom^2'), symmetry=1, barrier=(8.99978,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.10238,'amu*angstrom^2'), symmetry=1, barrier=(8.99978,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.10238,'amu*angstrom^2'), symmetry=1, barrier=(8.99978,'kJ/mol'), semiclassical=False),
    ],
    spinMultiplicity = 4,
    opticalIsomers = 1,
    molecularWeight = (117.08,'amu'),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[0.589122,0.0833615,-0.000133118,1.18835e-07,-4.19642e-11,16672.6,34.3041], Tmin=(100,'K'), Tmax=(815.942,'K')), NASAPolynomial(coeffs=[7.85573,0.0339726,-1.70166e-05,3.29823e-09,-2.29107e-13,15945,3.53298], Tmin=(815.942,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(137.67,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(291.007,'J/(mol*K)'), comment="""Thermo group additivity estimation: group(O2s-OsCs) + group(O2s-CsH) + group(O2s-CsH) + group(O2s-OsH) + group(Cs-(Cds-Cds)CsOsH) + group(Cs-CsOsOsH) + group(Cds-CdsCsH) + group(Cds-CdsHH) + radical(CCOJ) + radical(ROOJ) + radical(Cds_S)"""),
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
    label = 'C=C=C(O)C([O])O[O](28905)',
    structure = SMILES('C=C=C(O)C([O])O[O]'),
    E0 = (21.3002,'kJ/mol'),
    modes = [
        HarmonicOscillator(frequencies=([492.5,1135,1000,3615,1277.5,1000,540,610,2055,350,440,435,1725,2950,3100,1380,975,1025,1650,1380,1390,370,380,2900,435,180,180],'cm^-1')),
        HinderedRotor(inertia=(0.570147,'amu*angstrom^2'), symmetry=1, barrier=(13.1088,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.570112,'amu*angstrom^2'), symmetry=1, barrier=(13.108,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.570081,'amu*angstrom^2'), symmetry=1, barrier=(13.1073,'kJ/mol'), semiclassical=False),
    ],
    spinMultiplicity = 3,
    opticalIsomers = 1,
    molecularWeight = (116.072,'amu'),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[0.345631,0.0863753,-0.000135224,1.1129e-07,-3.61481e-11,2687.69,29.445], Tmin=(100,'K'), Tmax=(812.015,'K')), NASAPolynomial(coeffs=[11.7492,0.0249762,-1.21524e-05,2.32307e-09,-1.60086e-13,1007.99,-22.1352], Tmin=(812.015,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(21.3002,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(270.22,'J/(mol*K)'), comment="""Thermo group additivity estimation: group(O2s-OsCs) + group(O2s-CsH) + group(O2s-(Cds-Cd)H) + group(O2s-OsH) + group(Cs-CsOsOsH) + group(Cds-CdsCsOs) + group(Cds-CdsHH) + group(Cdd-CdsCds) + radical(ROOJ) + radical(CCOJ)"""),
)

species(
    label = 'C=[C]C(O)C(=O)O[O](28906)',
    structure = SMILES('C=[C]C(O)C(=O)O[O]'),
    E0 = (-105.954,'kJ/mol'),
    modes = [
        HarmonicOscillator(frequencies=([1685,370,492.5,1135,1000,2950,3100,1380,975,1025,1650,3615,1277.5,1000,1380,1390,370,380,2900,435,200,800,1000,1200,1400,1600],'cm^-1')),
        HinderedRotor(inertia=(0.156089,'amu*angstrom^2'), symmetry=1, barrier=(3.5888,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.156089,'amu*angstrom^2'), symmetry=1, barrier=(3.5888,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.156089,'amu*angstrom^2'), symmetry=1, barrier=(3.5888,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.156089,'amu*angstrom^2'), symmetry=1, barrier=(3.5888,'kJ/mol'), semiclassical=False),
    ],
    spinMultiplicity = 3,
    opticalIsomers = 1,
    molecularWeight = (116.072,'amu'),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[0.899874,0.0583968,-4.0534e-05,4.31238e-09,4.05739e-12,-12623.1,32.154], Tmin=(100,'K'), Tmax=(1004.72,'K')), NASAPolynomial(coeffs=[16.8218,0.0150497,-5.7399e-06,1.07902e-09,-7.82056e-14,-16834.1,-49.7691], Tmin=(1004.72,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(-105.954,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(266.063,'J/(mol*K)'), comment="""Thermo group additivity estimation: group(O2s-CsH) + group(O2s-O2s(Cds-O2d)) + group(O2s-OsH) + group(Cs-(Cds-O2d)CsOsH) + group(Cds-CdsCsH) + group(Cds-OdCsOs) + group(Cds-CdsHH) + radical(C(=O)OOJ) + radical(Cds_S)"""),
)

species(
    label = 'C#CC(O)C([O])O[O](28907)',
    structure = SMILES('C#CC(O)C([O])O[O]'),
    E0 = (65.7779,'kJ/mol'),
    modes = [
        HarmonicOscillator(frequencies=([492.5,1135,1000,750,770,3400,2100,3615,1277.5,1000,2175,525,1380,1383.33,1386.67,1390,370,373.333,376.667,380,2800,3000,430,440,222.506,3484.39],'cm^-1')),
        HinderedRotor(inertia=(0.307712,'amu*angstrom^2'), symmetry=1, barrier=(10.8066,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.307627,'amu*angstrom^2'), symmetry=1, barrier=(10.8069,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(1.37667,'amu*angstrom^2'), symmetry=1, barrier=(48.3693,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(1.37699,'amu*angstrom^2'), symmetry=1, barrier=(48.3694,'kJ/mol'), semiclassical=False),
    ],
    spinMultiplicity = 3,
    opticalIsomers = 1,
    molecularWeight = (116.072,'amu'),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[0.550155,0.0821488,-0.000128755,1.07583e-07,-3.51715e-11,8029.49,31.5892], Tmin=(100,'K'), Tmax=(852.933,'K')), NASAPolynomial(coeffs=[10.3122,0.026083,-1.2068e-05,2.24074e-09,-1.51245e-13,6738.33,-11.7608], Tmin=(852.933,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(65.7779,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(266.063,'J/(mol*K)'), comment="""Thermo group additivity estimation: group(O2s-OsCs) + group(O2s-CsH) + group(O2s-CsH) + group(O2s-OsH) + group(Cs-(Cds-Cds)CsOsH) + group(Cs-CsOsOsH) + group(Ct-CtCs) + group(Ct-CtH) + radical(ROOJ) + radical(CCOJ)"""),
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
    label = 'C=C=CC([O])O[O](18038)',
    structure = SMILES('C=C=CC([O])O[O]'),
    E0 = (235.859,'kJ/mol'),
    modes = [
        HarmonicOscillator(frequencies=([492.5,1135,1000,2950,3100,1380,975,1025,1650,1380,1390,370,380,2900,435,540,610,2055,3010,987.5,1337.5,450,1655,180,180],'cm^-1')),
        HinderedRotor(inertia=(0.362741,'amu*angstrom^2'), symmetry=1, barrier=(8.34012,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.362571,'amu*angstrom^2'), symmetry=1, barrier=(8.33621,'kJ/mol'), semiclassical=False),
    ],
    spinMultiplicity = 3,
    opticalIsomers = 1,
    molecularWeight = (100.073,'amu'),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[1.32638,0.0679542,-0.000114126,1.0911e-07,-4.04349e-11,28454.7,25.8633], Tmin=(100,'K'), Tmax=(831.164,'K')), NASAPolynomial(coeffs=[4.52737,0.0333134,-1.68949e-05,3.27702e-09,-2.2711e-13,28587,15.0094], Tmin=(831.164,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(235.859,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(249.434,'J/(mol*K)'), comment="""Thermo group additivity estimation: group(O2s-OsCs) + group(O2s-CsH) + group(O2s-OsH) + group(Cs-CsOsOsH) + group(Cds-CdsCsH) + group(Cds-CdsHH) + group(Cdd-CdsCds) + radical(CCOJ) + radical(ROOJ)"""),
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
    label = 'C=[C]C(O)C=O(27744)',
    structure = SMILES('C=[C]C(O)C=O'),
    E0 = (-32.4957,'kJ/mol'),
    modes = [
        HarmonicOscillator(frequencies=([1685,370,2782.5,750,1395,475,1775,1000,2950,3100,1380,975,1025,1650,3615,1277.5,1000,1380,1390,370,380,2900,435,616.43],'cm^-1')),
        HinderedRotor(inertia=(0.478592,'amu*angstrom^2'), symmetry=1, barrier=(11.0038,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.0407669,'amu*angstrom^2'), symmetry=1, barrier=(10.9798,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.0409213,'amu*angstrom^2'), symmetry=1, barrier=(10.9972,'kJ/mol'), semiclassical=False),
    ],
    spinMultiplicity = 2,
    opticalIsomers = 1,
    molecularWeight = (85.0813,'amu'),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[1.78798,0.043671,-3.37416e-05,1.31659e-08,-2.06728e-12,-3824.66,24.676], Tmin=(100,'K'), Tmax=(1501.82,'K')), NASAPolynomial(coeffs=[11.867,0.0168258,-6.92882e-06,1.26348e-09,-8.59241e-14,-6852.02,-28.0484], Tmin=(1501.82,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(-32.4957,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(245.277,'J/(mol*K)'), comment="""Thermo group additivity estimation: group(O2s-CsH) + group(Cs-(Cds-O2d)CsOsH) + group(Cds-CdsCsH) + group(Cds-OdCsH) + group(Cds-CdsHH) + radical(Cds_S)"""),
)

species(
    label = '[CH2]C=C(O)C([O])O[O](15274)',
    structure = SMILES('[CH2]C=C(O)C([O])O[O]'),
    E0 = (-3.80495,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
    molecularWeight = (117.08,'amu'),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[0.469743,0.0821335,-0.000109379,7.73689e-08,-2.19144e-11,-334.476,30.2376], Tmin=(100,'K'), Tmax=(861.94,'K')), NASAPolynomial(coeffs=[12.4709,0.0264372,-1.24491e-05,2.39514e-09,-1.67827e-13,-2403.25,-25.8777], Tmin=(861.94,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(-3.80495,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(291.007,'J/(mol*K)'), comment="""Thermo group additivity estimation: group(O2s-OsCs) + group(O2s-CsH) + group(O2s-(Cds-Cd)H) + group(O2s-OsH) + group(Cs-CsOsOsH) + group(Cs-(Cds-Cds)HHH) + group(Cds-CdsCsOs) + group(Cds-CdsCsH) + radical(CCOJ) + radical(Allyl_P) + radical(ROOJ)"""),
)

species(
    label = 'C=[C]C(O)[C](O)O[O](28908)',
    structure = SMILES('C=[C]C(O)[C](O)O[O]'),
    E0 = (117.211,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
    molecularWeight = (117.08,'amu'),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[0.342386,0.0858088,-0.000127078,9.95613e-08,-3.1017e-11,14223.9,35.0239], Tmin=(100,'K'), Tmax=(787.68,'K')), NASAPolynomial(coeffs=[12.2249,0.0254734,-1.2192e-05,2.33577e-09,-1.62048e-13,12351.8,-19.4677], Tmin=(787.68,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(117.211,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(286.849,'J/(mol*K)'), comment="""Thermo group additivity estimation: group(O2s-OsCs) + group(O2s-CsH) + group(O2s-CsH) + group(O2s-OsH) + group(Cs-(Cds-Cds)CsOsH) + group(Cs-CsOsOsH) + group(Cds-CdsCsH) + group(Cds-CdsHH) + radical(Cs_P) + radical(ROOJ) + radical(Cds_S)"""),
)

species(
    label = '[CH]=CC(O)C([O])O[O](15282)',
    structure = SMILES('[CH]=CC(O)C([O])O[O]'),
    E0 = (146.924,'kJ/mol'),
    modes = [
        HarmonicOscillator(frequencies=([3120,650,792.5,1650,492.5,1135,1000,3010,987.5,1337.5,450,1655,3615,1277.5,1000,1380,1383.33,1386.67,1390,370,373.333,376.667,380,2800,3000,430,440,265.773,265.854],'cm^-1')),
        HinderedRotor(inertia=(0.00238562,'amu*angstrom^2'), symmetry=1, barrier=(0.119627,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.141589,'amu*angstrom^2'), symmetry=1, barrier=(7.10093,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.141498,'amu*angstrom^2'), symmetry=1, barrier=(7.10152,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.141509,'amu*angstrom^2'), symmetry=1, barrier=(7.10083,'kJ/mol'), semiclassical=False),
    ],
    spinMultiplicity = 4,
    opticalIsomers = 1,
    molecularWeight = (117.08,'amu'),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[0.578764,0.0818912,-0.000123448,1.04727e-07,-3.58781e-11,17787.7,34.273], Tmin=(100,'K'), Tmax=(781.395,'K')), NASAPolynomial(coeffs=[9.19068,0.0318127,-1.58097e-05,3.07611e-09,-2.15096e-13,16624.8,-3.97934], Tmin=(781.395,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(146.924,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(291.007,'J/(mol*K)'), comment="""Thermo group additivity estimation: group(O2s-OsCs) + group(O2s-CsH) + group(O2s-CsH) + group(O2s-OsH) + group(Cs-(Cds-Cds)CsOsH) + group(Cs-CsOsOsH) + group(Cds-CdsCsH) + group(Cds-CdsHH) + radical(CCOJ) + radical(ROOJ) + radical(Cds_P)"""),
)

species(
    label = 'C=[C]C(O)[C]([O])OO(28909)',
    structure = SMILES('C=[C]C(O)[C]([O])OO'),
    E0 = (190.911,'kJ/mol'),
    modes = [
        HarmonicOscillator(frequencies=([360,370,350,3615,1277.5,1000,2950,3100,1380,975,1025,1650,1380,1390,370,380,2900,435,3615,1310,387.5,850,1000,1685,370,200,800,1600],'cm^-1')),
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
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[0.384061,0.0874302,-0.000138009,1.20866e-07,-4.24205e-11,23083.9,34.9695], Tmin=(100,'K'), Tmax=(783.615,'K')), NASAPolynomial(coeffs=[9.18777,0.0334294,-1.72938e-05,3.40965e-09,-2.39754e-13,21982.4,-3.58133], Tmin=(783.615,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(190.911,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(286.849,'J/(mol*K)'), comment="""Thermo group additivity estimation: group(O2s-OsCs) + group(O2s-CsH) + group(O2s-CsH) + group(O2s-OsH) + group(Cs-(Cds-Cds)CsOsH) + group(Cs-CsOsOsH) + group(Cds-CdsCsH) + group(Cds-CdsHH) + radical(Cs_P) + radical(CCOJ) + radical(Cds_S)"""),
)

species(
    label = 'C=CC([O])C([O])O[O](12832)',
    structure = SMILES('C=CC([O])C([O])O[O]'),
    E0 = (130.189,'kJ/mol'),
    modes = [
        HarmonicOscillator(frequencies=([492.5,1135,1000,2950,3100,1380,975,1025,1650,1380,1383.33,1386.67,1390,370,373.333,376.667,380,2800,3000,430,440,3010,987.5,1337.5,450,1655,303.693,303.698,303.973,1619.48],'cm^-1')),
        HinderedRotor(inertia=(0.0018236,'amu*angstrom^2'), symmetry=1, barrier=(0.119627,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.169329,'amu*angstrom^2'), symmetry=1, barrier=(11.1785,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.16935,'amu*angstrom^2'), symmetry=1, barrier=(11.1766,'kJ/mol'), semiclassical=False),
    ],
    spinMultiplicity = 4,
    opticalIsomers = 1,
    molecularWeight = (117.08,'amu'),
    collisionModel = TransportData(shapeIndex=2, epsilon=(4795.99,'J/mol'), sigma=(7.57616,'angstroms'), dipoleMoment=(0,'C*m'), polarizability=(0,'angstroms^3'), rotrelaxcollnum=0, comment="""Epsilon & sigma estimated with Tc=749.12 K, Pc=25.03 bar (from Joback method)"""),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[1.00385,0.0719129,-9.10648e-05,6.76519e-08,-2.11444e-11,15760.5,32.1676], Tmin=(100,'K'), Tmax=(769.121,'K')), NASAPolynomial(coeffs=[8.31821,0.0338708,-1.68683e-05,3.33572e-09,-2.37545e-13,14635.4,-1.19955], Tmin=(769.121,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(130.189,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(295.164,'J/(mol*K)'), comment="""Thermo group additivity estimation: group(O2s-OsCs) + group(O2s-CsH) + group(O2s-CsH) + group(O2s-OsH) + group(Cs-(Cds-Cds)CsOsH) + group(Cs-CsOsOsH) + group(Cds-CdsCsH) + group(Cds-CdsHH) + radical(ROOJ) + radical(CCOJ) + radical(CC(C)OJ)"""),
)

species(
    label = 'C=CC(O)[C]([O])O[O](15277)',
    structure = SMILES('C=CC(O)[C]([O])O[O]'),
    E0 = (105.074,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
    molecularWeight = (117.08,'amu'),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[0.704436,0.0789267,-0.000116634,9.88545e-08,-3.41263e-11,12750,34.107], Tmin=(100,'K'), Tmax=(769.753,'K')), NASAPolynomial(coeffs=[8.64153,0.0325293,-1.61795e-05,3.15715e-09,-2.21435e-13,11680.8,-1.11635], Tmin=(769.753,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(105.074,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(291.007,'J/(mol*K)'), comment="""Thermo group additivity estimation: group(O2s-OsCs) + group(O2s-CsH) + group(O2s-CsH) + group(O2s-OsH) + group(Cs-(Cds-Cds)CsOsH) + group(Cs-CsOsOsH) + group(Cds-CdsCsH) + group(Cds-CdsHH) + radical(ROOJ) + radical(CCOJ) + radical(Cs_P)"""),
)

species(
    label = '[CH2][C]=C(O)C(O)O[O](28910)',
    structure = SMILES('[CH2][C]=C(O)C(O)O[O]'),
    E0 = (8.33168,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
    molecularWeight = (117.08,'amu'),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[-0.0148528,0.0906326,-0.000126384,8.81151e-08,-2.3903e-11,1144.61,31.5836], Tmin=(100,'K'), Tmax=(909.689,'K')), NASAPolynomial(coeffs=[16.2845,0.0189638,-8.20956e-06,1.51224e-09,-1.03211e-13,-1820.9,-45.5087], Tmin=(909.689,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(8.33168,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(286.849,'J/(mol*K)'), comment="""Thermo group additivity estimation: group(O2s-OsCs) + group(O2s-CsH) + group(O2s-(Cds-Cd)H) + group(O2s-OsH) + group(Cs-CsOsOsH) + group(Cs-(Cds-Cds)HHH) + group(Cds-CdsCsOs) + group(Cds-CdsCsH) + radical(Allyl_P) + radical(Cds_S) + radical(ROOJ)"""),
)

species(
    label = 'C=[C]C([O])C(O)O[O](15280)',
    structure = SMILES('C=[C]C([O])C(O)O[O]'),
    E0 = (142.325,'kJ/mol'),
    modes = [
        HarmonicOscillator(frequencies=([1685,370,492.5,1135,1000,2950,3100,1380,975,1025,1650,3615,1277.5,1000,1380,1383.33,1386.67,1390,370,373.333,376.667,380,2800,3000,430,440,328.989,329.031,329.106],'cm^-1')),
        HinderedRotor(inertia=(0.00155636,'amu*angstrom^2'), symmetry=1, barrier=(0.119627,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.136515,'amu*angstrom^2'), symmetry=1, barrier=(10.4923,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.136698,'amu*angstrom^2'), symmetry=1, barrier=(10.4913,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.136587,'amu*angstrom^2'), symmetry=1, barrier=(10.4917,'kJ/mol'), semiclassical=False),
    ],
    spinMultiplicity = 4,
    opticalIsomers = 1,
    molecularWeight = (117.08,'amu'),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[0.603279,0.0793657,-0.000104119,7.27578e-08,-2.04397e-11,17236,33.2157], Tmin=(100,'K'), Tmax=(867.107,'K')), NASAPolynomial(coeffs=[12.0597,0.0265166,-1.26959e-05,2.46833e-09,-1.74191e-13,15249.2,-20.4215], Tmin=(867.107,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(142.325,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(291.007,'J/(mol*K)'), comment="""Thermo group additivity estimation: group(O2s-OsCs) + group(O2s-CsH) + group(O2s-CsH) + group(O2s-OsH) + group(Cs-(Cds-Cds)CsOsH) + group(Cs-CsOsOsH) + group(Cds-CdsCsH) + group(Cds-CdsHH) + radical(ROOJ) + radical(Cds_S) + radical(CC(C)OJ)"""),
)

species(
    label = '[CH2][C]=C(O)C([O])OO(28911)',
    structure = SMILES('[CH2][C]=C(O)C([O])OO'),
    E0 = (82.0321,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
    molecularWeight = (117.08,'amu'),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[0.156036,0.0905254,-0.000130199,9.84138e-08,-2.9682e-11,9999.16,31.0783], Tmin=(100,'K'), Tmax=(811.564,'K')), NASAPolynomial(coeffs=[12.9514,0.027461,-1.364e-05,2.66668e-09,-1.87786e-13,7922.28,-27.9805], Tmin=(811.564,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(82.0321,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(286.849,'J/(mol*K)'), comment="""Thermo group additivity estimation: group(O2s-OsCs) + group(O2s-CsH) + group(O2s-(Cds-Cd)H) + group(O2s-OsH) + group(Cs-CsOsOsH) + group(Cs-(Cds-Cds)HHH) + group(Cds-CdsCsOs) + group(Cds-CdsCsH) + radical(CCOJ) + radical(Allyl_P) + radical(Cds_S)"""),
)

species(
    label = '[CH]=[C]C(O)C(O)O[O](28912)',
    structure = SMILES('[CH]=[C]C(O)C(O)O[O]'),
    E0 = (159.061,'kJ/mol'),
    modes = [
        HarmonicOscillator(frequencies=([1685,370,492.5,1135,1000,3580,3650,1210,1345,900,1100,3120,650,792.5,1650,1380,1383.33,1386.67,1390,370,373.333,376.667,380,2800,3000,430,440,180],'cm^-1')),
        HinderedRotor(inertia=(0.545629,'amu*angstrom^2'), symmetry=1, barrier=(12.5451,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.00110466,'amu*angstrom^2'), symmetry=1, barrier=(12.5424,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.545523,'amu*angstrom^2'), symmetry=1, barrier=(12.5426,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.974759,'amu*angstrom^2'), symmetry=1, barrier=(22.4116,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.545632,'amu*angstrom^2'), symmetry=1, barrier=(12.5452,'kJ/mol'), semiclassical=False),
    ],
    spinMultiplicity = 4,
    opticalIsomers = 1,
    molecularWeight = (117.08,'amu'),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[0.335461,0.087101,-0.000126478,9.28243e-08,-2.55768e-11,19256.5,34.7796], Tmin=(100,'K'), Tmax=(685.501,'K')), NASAPolynomial(coeffs=[12.6611,0.0249716,-1.19561e-05,2.28815e-09,-1.58596e-13,17336.6,-21.7087], Tmin=(685.501,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(159.061,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(286.849,'J/(mol*K)'), comment="""Thermo group additivity estimation: group(O2s-OsCs) + group(O2s-CsH) + group(O2s-CsH) + group(O2s-OsH) + group(Cs-(Cds-Cds)CsOsH) + group(Cs-CsOsOsH) + group(Cds-CdsCsH) + group(Cds-CdsHH) + radical(ROOJ) + radical(Cds_S) + radical(Cds_P)"""),
)

species(
    label = 'C=[C]C([O])C([O])OO(15284)',
    structure = SMILES('C=[C]C([O])C([O])OO'),
    E0 = (216.026,'kJ/mol'),
    modes = [
        HarmonicOscillator(frequencies=([1685,370,3615,1310,387.5,850,1000,1380,1383.33,1386.67,1390,370,373.333,376.667,380,2800,3000,430,440,2950,3100,1380,975,1025,1650,386.396,386.396,386.398,1857.56],'cm^-1')),
        HinderedRotor(inertia=(0.0753911,'amu*angstrom^2'), symmetry=1, barrier=(7.98754,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.0753911,'amu*angstrom^2'), symmetry=1, barrier=(7.98754,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.2713,'amu*angstrom^2'), symmetry=1, barrier=(28.7436,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.271299,'amu*angstrom^2'), symmetry=1, barrier=(28.7436,'kJ/mol'), semiclassical=False),
    ],
    spinMultiplicity = 4,
    opticalIsomers = 1,
    molecularWeight = (117.08,'amu'),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[0.833509,0.0780859,-0.000100765,6.69715e-08,-1.45836e-11,26088,32.5223], Tmin=(100,'K'), Tmax=(595.316,'K')), NASAPolynomial(coeffs=[8.80344,0.0349002,-1.80687e-05,3.6107e-09,-2.57867e-13,24955.4,-3.33693], Tmin=(595.316,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(216.026,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(291.007,'J/(mol*K)'), comment="""Thermo group additivity estimation: group(O2s-OsCs) + group(O2s-CsH) + group(O2s-CsH) + group(O2s-OsH) + group(Cs-(Cds-Cds)CsOsH) + group(Cs-CsOsOsH) + group(Cds-CdsCsH) + group(Cds-CdsHH) + radical(CCOJ) + radical(CC(C)OJ) + radical(Cds_S)"""),
)

species(
    label = '[CH]=[C]C(O)C([O])OO(28913)',
    structure = SMILES('[CH]=[C]C(O)C([O])OO'),
    E0 = (232.761,'kJ/mol'),
    modes = [
        HarmonicOscillator(frequencies=([1685,370,3615,1310,387.5,850,1000,3615,1277.5,1000,3120,650,792.5,1650,1380,1383.33,1386.67,1390,370,373.333,376.667,380,2800,3000,430,440,205.892,1827.45],'cm^-1')),
        HinderedRotor(inertia=(0.333271,'amu*angstrom^2'), symmetry=1, barrier=(10.0249,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.333255,'amu*angstrom^2'), symmetry=1, barrier=(10.0249,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.906688,'amu*angstrom^2'), symmetry=1, barrier=(27.2746,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.00397645,'amu*angstrom^2'), symmetry=1, barrier=(0.119627,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.906702,'amu*angstrom^2'), symmetry=1, barrier=(27.2746,'kJ/mol'), semiclassical=False),
    ],
    spinMultiplicity = 4,
    opticalIsomers = 1,
    molecularWeight = (117.08,'amu'),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[0.261927,0.0903477,-0.000144632,1.26449e-07,-4.4028e-11,28121.5,35.1232], Tmin=(100,'K'), Tmax=(792.127,'K')), NASAPolynomial(coeffs=[9.72837,0.0327285,-1.69336e-05,3.33096e-09,-2.33616e-13,26929.7,-6.39692], Tmin=(792.127,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(232.761,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(286.849,'J/(mol*K)'), comment="""Thermo group additivity estimation: group(O2s-OsCs) + group(O2s-CsH) + group(O2s-CsH) + group(O2s-OsH) + group(Cs-(Cds-Cds)CsOsH) + group(Cs-CsOsOsH) + group(Cds-CdsCsH) + group(Cds-CdsHH) + radical(Cds_P) + radical(Cds_S) + radical(CCOJ)"""),
)

species(
    label = 'C=[C][CH]C([O])O[O](18045)',
    structure = SMILES('C=[C][CH]C([O])O[O]'),
    E0 = (427.393,'kJ/mol'),
    modes = [
        HarmonicOscillator(frequencies=([1685,370,492.5,1135,1000,2950,3100,1380,975,1025,1650,1380,1390,370,380,2900,435,3025,407.5,1350,352.5,180,378.242,4000],'cm^-1')),
        HinderedRotor(inertia=(1.44583,'amu*angstrom^2'), symmetry=1, barrier=(33.2424,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(2.03979,'amu*angstrom^2'), symmetry=1, barrier=(46.8989,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(2.04739,'amu*angstrom^2'), symmetry=1, barrier=(47.0735,'kJ/mol'), semiclassical=False),
    ],
    spinMultiplicity = 5,
    opticalIsomers = 1,
    molecularWeight = (100.073,'amu'),
    collisionModel = TransportData(shapeIndex=2, epsilon=(4216.54,'J/mol'), sigma=(6.85526,'angstroms'), dipoleMoment=(0,'C*m'), polarizability=(0,'angstroms^3'), rotrelaxcollnum=0, comment="""Epsilon & sigma estimated with Tc=658.61 K, Pc=29.7 bar (from Joback method)"""),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[1.05198,0.0729387,-0.000118694,1.08795e-07,-3.93024e-11,51501.8,26.5742], Tmin=(100,'K'), Tmax=(815.932,'K')), NASAPolynomial(coeffs=[6.46296,0.0318145,-1.62559e-05,3.16944e-09,-2.20735e-13,51104.7,4.54774], Tmin=(815.932,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(427.393,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(245.277,'J/(mol*K)'), comment="""Thermo group additivity estimation: group(O2s-OsCs) + group(O2s-CsH) + group(O2s-OsH) + group(Cs-(Cds-Cds)CsHH) + group(Cs-CsOsOsH) + group(Cds-CdsCsH) + group(Cds-CdsHH) + radical(ROOJ) + radical(CCOJ) + radical(Cds_S) + radical(C=CCJCO)"""),
)

species(
    label = 'C=[C]C([O])C([O])O[O](15288)',
    structure = SMILES('C=[C]C([O])C([O])O[O]'),
    E0 = (368.03,'kJ/mol'),
    modes = [
        HarmonicOscillator(frequencies=([1685,370,492.5,1135,1000,1380,1383.33,1386.67,1390,370,373.333,376.667,380,2800,3000,430,440,2950,3100,1380,975,1025,1650,370.51,370.965,371.043,1581.14],'cm^-1')),
        HinderedRotor(inertia=(0.082548,'amu*angstrom^2'), symmetry=1, barrier=(8.00413,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.0831752,'amu*angstrom^2'), symmetry=1, barrier=(8.03029,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.0819688,'amu*angstrom^2'), symmetry=1, barrier=(8.03395,'kJ/mol'), semiclassical=False),
    ],
    spinMultiplicity = 5,
    opticalIsomers = 1,
    molecularWeight = (116.072,'amu'),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[0.752342,0.0793718,-0.000127287,1.14031e-07,-4.05434e-11,44373.1,33.5259], Tmin=(100,'K'), Tmax=(805.589,'K')), NASAPolynomial(coeffs=[7.74161,0.0322653,-1.64811e-05,3.22281e-09,-2.25133e-13,43649.5,3.81553], Tmin=(805.589,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(368.03,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(270.22,'J/(mol*K)'), comment="""Thermo group additivity estimation: group(O2s-OsCs) + group(O2s-CsH) + group(O2s-CsH) + group(O2s-OsH) + group(Cs-(Cds-Cds)CsOsH) + group(Cs-CsOsOsH) + group(Cds-CdsCsH) + group(Cds-CdsHH) + radical(ROOJ) + radical(Cds_S) + radical(CC(C)OJ) + radical(CCOJ)"""),
)

species(
    label = 'C=[C]C(O)[CH][O](28102)',
    structure = SMILES('C=[C]C(O)[CH][O]'),
    E0 = (294.689,'kJ/mol'),
    modes = [
        HarmonicOscillator(frequencies=([1685,370,3025,407.5,1350,352.5,2950,3100,1380,975,1025,1650,3615,1277.5,1000,1380,1390,370,380,2900,435,276.259,279.06,2001.85],'cm^-1')),
        HinderedRotor(inertia=(0.00227246,'amu*angstrom^2'), symmetry=1, barrier=(0.119627,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.191772,'amu*angstrom^2'), symmetry=1, barrier=(10.9158,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.189866,'amu*angstrom^2'), symmetry=1, barrier=(10.8672,'kJ/mol'), semiclassical=False),
    ],
    spinMultiplicity = 4,
    opticalIsomers = 1,
    molecularWeight = (85.0813,'amu'),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[1.40118,0.0638678,-0.000101269,9.07849e-08,-3.19727e-11,35530,26.2595], Tmin=(100,'K'), Tmax=(835.883,'K')), NASAPolynomial(coeffs=[6.37686,0.0274382,-1.32505e-05,2.52384e-09,-1.73332e-13,35139,5.78378], Tmin=(835.883,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(294.689,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(245.277,'J/(mol*K)'), comment="""Thermo group additivity estimation: group(O2s-CsH) + group(O2s-CsH) + group(Cs-(Cds-Cds)CsOsH) + group(Cs-CsOsHH) + group(Cds-CdsCsH) + group(Cds-CdsHH) + radical(CCOJ) + radical(Cds_S) + radical(CCsJOH)"""),
)

species(
    label = '[CH2][C]=C(O)C([O])O[O](28914)',
    structure = SMILES('[CH2][C]=C(O)C([O])O[O]'),
    E0 = (234.037,'kJ/mol'),
    modes = [
        HarmonicOscillator(frequencies=([492.5,1135,1000,3615,1277.5,1000,350,440,435,1725,3000,3100,440,815,1455,1000,1685,370,1380,1390,370,380,2900,435,216.811,2426.84],'cm^-1')),
        HinderedRotor(inertia=(0.372556,'amu*angstrom^2'), symmetry=1, barrier=(12.427,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.372547,'amu*angstrom^2'), symmetry=1, barrier=(12.427,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(1.02959,'amu*angstrom^2'), symmetry=1, barrier=(34.3439,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.372544,'amu*angstrom^2'), symmetry=1, barrier=(12.427,'kJ/mol'), semiclassical=False),
    ],
    spinMultiplicity = 5,
    opticalIsomers = 1,
    molecularWeight = (116.072,'amu'),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[0.254978,0.089061,-0.000143234,1.19865e-07,-3.92491e-11,28276.6,31.4702], Tmin=(100,'K'), Tmax=(832.869,'K')), NASAPolynomial(coeffs=[11.733,0.0251336,-1.2248e-05,2.32839e-09,-1.59384e-13,26670,-19.9726], Tmin=(832.869,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(234.037,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(266.063,'J/(mol*K)'), comment="""Thermo group additivity estimation: group(O2s-OsCs) + group(O2s-CsH) + group(O2s-(Cds-Cd)H) + group(O2s-OsH) + group(Cs-CsOsOsH) + group(Cs-(Cds-Cds)HHH) + group(Cds-CdsCsOs) + group(Cds-CdsCsH) + radical(CCOJ) + radical(ROOJ) + radical(Cds_S) + radical(Allyl_P)"""),
)

species(
    label = 'C=[C]C(O)[C]([O])O[O](28915)',
    structure = SMILES('C=[C]C(O)[C]([O])O[O]'),
    E0 = (342.916,'kJ/mol'),
    modes = [
        HarmonicOscillator(frequencies=([492.5,1135,1000,3615,1277.5,1000,360,370,350,2950,3100,1380,975,1025,1650,1380,1390,370,380,2900,435,1685,370,289.335,289.532,1969.7],'cm^-1')),
        HinderedRotor(inertia=(0.00201333,'amu*angstrom^2'), symmetry=1, barrier=(0.119627,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.137011,'amu*angstrom^2'), symmetry=1, barrier=(8.14948,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.137104,'amu*angstrom^2'), symmetry=1, barrier=(8.15083,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.137155,'amu*angstrom^2'), symmetry=1, barrier=(8.15083,'kJ/mol'), semiclassical=False),
    ],
    spinMultiplicity = 5,
    opticalIsomers = 1,
    molecularWeight = (116.072,'amu'),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[0.605662,0.0843428,-0.000144432,1.32143e-07,-4.67899e-11,41356.2,34.9321], Tmin=(100,'K'), Tmax=(835.631,'K')), NASAPolynomial(coeffs=[7.75914,0.0314852,-1.61339e-05,3.12813e-09,-2.16185e-13,40810.6,5.59415], Tmin=(835.631,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(342.916,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(266.063,'J/(mol*K)'), comment="""Thermo group additivity estimation: group(O2s-OsCs) + group(O2s-CsH) + group(O2s-CsH) + group(O2s-OsH) + group(Cs-(Cds-Cds)CsOsH) + group(Cs-CsOsOsH) + group(Cds-CdsCsH) + group(Cds-CdsHH) + radical(Cds_S) + radical(Cs_P) + radical(CCOJ) + radical(ROOJ)"""),
)

species(
    label = '[CH]=[C]C(O)C([O])O[O](28916)',
    structure = SMILES('[CH]=[C]C(O)C([O])O[O]'),
    E0 = (384.766,'kJ/mol'),
    modes = [
        HarmonicOscillator(frequencies=([1685,370,492.5,1135,1000,3615,1277.5,1000,3120,650,792.5,1650,1380,1383.33,1386.67,1390,370,373.333,376.667,380,2800,3000,430,440,181.156,183.97],'cm^-1')),
        HinderedRotor(inertia=(0.193734,'amu*angstrom^2'), symmetry=1, barrier=(4.54148,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.191453,'amu*angstrom^2'), symmetry=1, barrier=(4.55348,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.192907,'amu*angstrom^2'), symmetry=1, barrier=(4.54823,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.192979,'amu*angstrom^2'), symmetry=1, barrier=(4.54731,'kJ/mol'), semiclassical=False),
    ],
    spinMultiplicity = 5,
    opticalIsomers = 1,
    molecularWeight = (116.072,'amu'),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[0.488061,0.0872047,-0.000150856,1.37477e-07,-4.8306e-11,46393.5,35.0697], Tmin=(100,'K'), Tmax=(840.828,'K')), NASAPolynomial(coeffs=[8.27629,0.0308263,-1.57988e-05,3.05555e-09,-2.10564e-13,45767,2.90915], Tmin=(840.828,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(384.766,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(266.063,'J/(mol*K)'), comment="""Thermo group additivity estimation: group(O2s-OsCs) + group(O2s-CsH) + group(O2s-CsH) + group(O2s-OsH) + group(Cs-(Cds-Cds)CsOsH) + group(Cs-CsOsOsH) + group(Cds-CdsCsH) + group(Cds-CdsHH) + radical(ROOJ) + radical(Cds_P) + radical(Cds_S) + radical(CCOJ)"""),
)

species(
    label = 'C=C=C(O)C(O)O[O](28917)',
    structure = SMILES('C=C=C(O)C(O)O[O]'),
    E0 = (-204.405,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
    molecularWeight = (117.08,'amu'),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[0.072034,0.0880022,-0.000118608,7.98722e-08,-2.09401e-11,-24444.1,29.5712], Tmin=(100,'K'), Tmax=(940.394,'K')), NASAPolynomial(coeffs=[16.3689,0.0186833,-8.04001e-06,1.48892e-09,-1.02385e-13,-27509.3,-48.0506], Tmin=(940.394,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(-204.405,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(291.007,'J/(mol*K)'), comment="""Thermo group additivity estimation: group(O2s-OsCs) + group(O2s-CsH) + group(O2s-(Cds-Cd)H) + group(O2s-OsH) + group(Cs-CsOsOsH) + group(Cds-CdsCsOs) + group(Cds-CdsHH) + group(Cdd-CdsCds) + radical(ROOJ)"""),
)

species(
    label = 'C=[C]C(O)C(=O)OO(27760)',
    structure = SMILES('C=[C]C(O)C(=O)OO'),
    E0 = (-300.343,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
    molecularWeight = (117.08,'amu'),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[0.61542,0.0622147,-3.65299e-05,-4.23949e-09,7.56954e-12,-35990.3,32.4168], Tmin=(100,'K'), Tmax=(1007.74,'K')), NASAPolynomial(coeffs=[18.7026,0.0162114,-6.44255e-06,1.25147e-09,-9.26609e-14,-40945.3,-61.4803], Tmin=(1007.74,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(-300.343,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(286.849,'J/(mol*K)'), comment="""Thermo group additivity estimation: group(O2s-CsH) + group(O2s-O2s(Cds-O2d)) + group(O2s-OsH) + group(Cs-(Cds-O2d)CsOsH) + group(Cds-CdsCsH) + group(Cds-OdCsOs) + group(Cds-CdsHH) + radical(Cds_S)"""),
)

species(
    label = 'C=CC(O)C(=O)O[O](15293)',
    structure = SMILES('C=CC(O)C(=O)O[O]'),
    E0 = (-343.796,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
    molecularWeight = (117.08,'amu'),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[0.944669,0.0536625,-1.52999e-05,-2.54941e-08,1.52312e-11,-41227,31.5192], Tmin=(100,'K'), Tmax=(975.653,'K')), NASAPolynomial(coeffs=[17.9073,0.0157235,-5.56127e-06,1.05319e-09,-7.87879e-14,-46041.1,-57.6063], Tmin=(975.653,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(-343.796,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(291.007,'J/(mol*K)'), comment="""Thermo group additivity estimation: group(O2s-CsH) + group(O2s-O2s(Cds-O2d)) + group(O2s-OsH) + group(Cs-(Cds-O2d)CsOsH) + group(Cds-CdsCsH) + group(Cds-OdCsOs) + group(Cds-CdsHH) + radical(C(=O)OOJ)"""),
)

species(
    label = 'C=C=C(O)C([O])OO(28918)',
    structure = SMILES('C=C=C(O)C([O])OO'),
    E0 = (-130.705,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
    molecularWeight = (117.08,'amu'),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[0.285619,0.0873452,-0.000120271,8.70208e-08,-2.51979e-11,-15591.4,28.9157], Tmin=(100,'K'), Tmax=(843.245,'K')), NASAPolynomial(coeffs=[12.928,0.0273726,-1.35847e-05,2.67097e-09,-1.8929e-13,-17723.4,-29.9203], Tmin=(843.245,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(-130.705,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(291.007,'J/(mol*K)'), comment="""Thermo group additivity estimation: group(O2s-OsCs) + group(O2s-CsH) + group(O2s-(Cds-Cd)H) + group(O2s-OsH) + group(Cs-CsOsOsH) + group(Cds-CdsCsOs) + group(Cds-CdsHH) + group(Cdd-CdsCds) + radical(CCOJ)"""),
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
    label = 'C=[C]C(O)C1OO1(28919)',
    structure = SMILES('C=[C]C(O)C1OO1'),
    E0 = (82.6658,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
    molecularWeight = (101.081,'amu'),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[1.24014,0.0481134,-1.09114e-05,-3.15287e-08,1.90241e-11,10053.4,27.1366], Tmin=(100,'K'), Tmax=(922.07,'K')), NASAPolynomial(coeffs=[17.6787,0.00955377,-1.46345e-06,1.62107e-10,-1.25792e-14,5629.64,-58.386], Tmin=(922.07,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(82.6658,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(270.22,'J/(mol*K)'), comment="""Thermo group additivity estimation: group(O2s-OsCs) + group(O2s-OsCs) + group(O2s-CsH) + group(Cs-(Cds-Cds)CsOsH) + group(Cs-CsOsOsH) + group(Cds-CdsCsH) + group(Cds-CdsHH) + ring(dioxirane) + radical(Cds_S)"""),
)

species(
    label = 'C=C1OC([O])C1O(28920)',
    structure = SMILES('C=C1OC([O])C1O'),
    E0 = (-202.57,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
    molecularWeight = (101.081,'amu'),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[1.46627,0.0402094,1.38903e-05,-5.84442e-08,2.89583e-11,-24257.9,20.6046], Tmin=(100,'K'), Tmax=(918.555,'K')), NASAPolynomial(coeffs=[18.0457,0.00850445,-4.59984e-07,-3.72963e-11,1.44755e-16,-29012,-67.2716], Tmin=(918.555,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(-202.57,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(278.535,'J/(mol*K)'), comment="""Thermo group additivity estimation: group(O2s-Cs(Cds-Cd)) + group(O2s-CsH) + group(O2s-CsH) + group(Cs-(Cds-Cds)CsOsH) + group(Cs-CsOsOsH) + group(Cds-CdsCsOs) + group(Cds-CdsHH) + ring(2methyleneoxetane) + radical(CCOJ)"""),
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
    label = 'C=[C]C(O)=C[O](28921)',
    structure = SMILES('C=[C]C(O)=C[O]'),
    E0 = (11.6938,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
    molecularWeight = (84.0734,'amu'),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[-0.226908,0.0685991,-7.82044e-05,4.03939e-08,-7.48367e-12,1579.65,23.2831], Tmin=(100,'K'), Tmax=(1612.29,'K')), NASAPolynomial(coeffs=[20.065,-0.000825657,4.13795e-06,-9.94445e-10,7.21889e-14,-2483.51,-76.6147], Tmin=(1612.29,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(11.6938,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(224.491,'J/(mol*K)'), comment="""Thermo group additivity estimation: group(O2s-(Cds-Cd)H) + group(O2s-(Cds-Cd)H) + group(Cds-Cds(Cds-Cds)O2s) + group(Cds-Cds(Cds-Cds)H) + group(Cds-CdsOsH) + group(Cds-CdsHH) + radical(C=CJC=C) + radical(C=COJ)"""),
)

species(
    label = 'C=[C]C(O)C1OOO1(28922)',
    structure = SMILES('C=[C]C(O)C1OOO1'),
    E0 = (123.878,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
    molecularWeight = (117.08,'amu'),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[1.14134,0.0441813,1.74963e-05,-6.09448e-08,2.75141e-11,15018.7,30.5876], Tmin=(100,'K'), Tmax=(983.373,'K')), NASAPolynomial(coeffs=[19.5235,0.0143559,-5.56895e-06,1.17136e-09,-9.38301e-14,9230.22,-68.8374], Tmin=(983.373,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(123.878,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(295.164,'J/(mol*K)'), comment="""Thermo group additivity estimation: group(O2s-OsCs) + group(O2s-OsCs) + group(O2s-CsH) + group(O2s-OsOs) + group(Cs-(Cds-Cds)CsOsH) + group(Cs-CsOsOsH) + group(Cds-CdsCsH) + group(Cds-CdsHH) + ring(Cyclobutane) + radical(Cds_S)"""),
)

species(
    label = 'C=C1OC(O[O])C1O(28884)',
    structure = SMILES('C=C1OC(O[O])C1O'),
    E0 = (-204.765,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
    molecularWeight = (117.08,'amu'),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[0.695922,0.0562455,-9.49901e-06,-4.60638e-08,2.73901e-11,-24493.2,24.0303], Tmin=(100,'K'), Tmax=(907.377,'K')), NASAPolynomial(coeffs=[22.2597,0.0053905,1.49447e-06,-4.50845e-10,3.02087e-14,-30226.3,-87.9341], Tmin=(907.377,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(-204.765,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(299.321,'J/(mol*K)'), comment="""Thermo group additivity estimation: group(O2s-Cs(Cds-Cd)) + group(O2s-OsCs) + group(O2s-CsH) + group(O2s-OsH) + group(Cs-(Cds-Cds)CsOsH) + group(Cs-CsOsOsH) + group(Cds-CdsCsOs) + group(Cds-CdsHH) + ring(2methyleneoxetane) + radical(ROOJ)"""),
)

species(
    label = 'C=C1OOC([O])C1O(28923)',
    structure = SMILES('C=C1OOC([O])C1O'),
    E0 = (-167.571,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
    molecularWeight = (117.08,'amu'),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[1.61099,0.0393308,1.06778e-05,-3.96313e-08,1.66742e-11,-20056.7,28.0155], Tmin=(100,'K'), Tmax=(1039.04,'K')), NASAPolynomial(coeffs=[13.6937,0.0233185,-1.02407e-05,2.04354e-09,-1.50906e-13,-24214.1,-38.6625], Tmin=(1039.04,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(-167.571,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(303.478,'J/(mol*K)'), comment="""Thermo group additivity estimation: group(O2s-OsCs) + group(O2s-O2s(Cds-Cd)) + group(O2s-CsH) + group(O2s-CsH) + group(Cs-(Cds-Cds)CsOsH) + group(Cs-CsOsOsH) + group(Cds-CdsCsOs) + group(Cds-CdsHH) + ring(Cyclopentane) + radical(CCOJ)"""),
)

species(
    label = 'C=[C][CH]C(O[O])OO(28924)',
    structure = SMILES('C=[C][CH]C(O[O])OO'),
    E0 = (273.193,'kJ/mol'),
    modes = [
        HarmonicOscillator(frequencies=([492.5,1135,1000,3025,407.5,1350,352.5,2950,3100,1380,975,1025,1650,1380,1390,370,380,2900,435,1685,370,3615,1310,387.5,850,1000,241.273,241.282],'cm^-1')),
        HinderedRotor(inertia=(0.383008,'amu*angstrom^2'), symmetry=1, barrier=(15.824,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(1.01494,'amu*angstrom^2'), symmetry=1, barrier=(41.9343,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.38307,'amu*angstrom^2'), symmetry=1, barrier=(15.824,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(1.01507,'amu*angstrom^2'), symmetry=1, barrier=(41.9345,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(1.01496,'amu*angstrom^2'), symmetry=1, barrier=(41.9343,'kJ/mol'), semiclassical=False),
    ],
    spinMultiplicity = 4,
    opticalIsomers = 1,
    molecularWeight = (117.08,'amu'),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[0.203392,0.090135,-0.000127655,9.73319e-08,-2.99627e-11,32988.2,29.5373], Tmin=(100,'K'), Tmax=(792.14,'K')), NASAPolynomial(coeffs=[11.9118,0.0310087,-1.56859e-05,3.09309e-09,-2.19025e-13,31133.4,-24.2199], Tmin=(792.14,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(273.193,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(286.849,'J/(mol*K)'), comment="""Thermo group additivity estimation: group(O2s-OsCs) + group(O2s-OsCs) + group(O2s-OsH) + group(O2s-OsH) + group(Cs-CsOsOsH) + group(Cs-(Cds-Cds)CsHH) + group(Cds-CdsCsH) + group(Cds-CdsHH) + radical(ROOJ) + radical(Cds_S) + radical(C=CCJCO)"""),
)

species(
    label = 'C=[C]C(O)C([O])[O](28925)',
    structure = SMILES('C=[C]C(O)C([O])[O]'),
    E0 = (139.865,'kJ/mol'),
    modes = [
        HarmonicOscillator(frequencies=([1685,370,2950,3100,1380,975,1025,1650,3615,1277.5,1000,1380,1383.33,1386.67,1390,370,373.333,376.667,380,2800,3000,430,440,264.106,264.566,1758.87,1759.36],'cm^-1')),
        HinderedRotor(inertia=(0.0024102,'amu*angstrom^2'), symmetry=1, barrier=(0.119627,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.145864,'amu*angstrom^2'), symmetry=1, barrier=(7.24453,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.145965,'amu*angstrom^2'), symmetry=1, barrier=(7.24474,'kJ/mol'), semiclassical=False),
    ],
    spinMultiplicity = 4,
    opticalIsomers = 1,
    molecularWeight = (101.081,'amu'),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[1.35663,0.0673854,-0.000110094,1.0723e-07,-4.09039e-11,16908.1,30.194], Tmin=(100,'K'), Tmax=(815.36,'K')), NASAPolynomial(coeffs=[3.57115,0.0372061,-1.90398e-05,3.72799e-09,-2.60513e-13,17189,23.8994], Tmin=(815.36,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(139.865,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(270.22,'J/(mol*K)'), comment="""Thermo group additivity estimation: group(O2s-CsH) + group(O2s-CsH) + group(O2s-CsH) + group(Cs-(Cds-Cds)CsOsH) + group(Cs-CsOsOsH) + group(Cds-CdsCsH) + group(Cds-CdsHH) + radical(CCOJ) + radical(CCOJ) + radical(Cds_S)"""),
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
    label = '[O]OC([O])[CH]O(8411)',
    structure = SMILES('[O]OC([O])[CH]O'),
    E0 = (12.8832,'kJ/mol'),
    modes = [
        HarmonicOscillator(frequencies=([3025,407.5,1350,352.5,492.5,1135,1000,3615,1277.5,1000,1380,1390,370,380,2900,435,180,1970.65],'cm^-1')),
        HinderedRotor(inertia=(0.295327,'amu*angstrom^2'), symmetry=1, barrier=(6.79016,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.295363,'amu*angstrom^2'), symmetry=1, barrier=(6.79097,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.295459,'amu*angstrom^2'), symmetry=1, barrier=(6.7932,'kJ/mol'), semiclassical=False),
    ],
    spinMultiplicity = 4,
    opticalIsomers = 1,
    molecularWeight = (91.0428,'amu'),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[1.40254,0.0697427,-0.000139515,1.37248e-07,-4.96232e-11,1631,24.9335], Tmin=(100,'K'), Tmax=(877.976,'K')), NASAPolynomial(coeffs=[4.40672,0.0260101,-1.34669e-05,2.55852e-09,-1.72166e-13,2261.51,17.4259], Tmin=(877.976,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(12.8832,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(195.39,'J/(mol*K)'), comment="""Thermo group additivity estimation: group(O2s-OsCs) + group(O2s-CsH) + group(O2s-CsH) + group(O2s-OsH) + group(Cs-CsOsOsH) + group(Cs-CsOsHH) + radical(CCsJOH) + radical(CCOJ) + radical(ROOJ)"""),
)

species(
    label = 'C=[C]C(O)[CH]O[O](28926)',
    structure = SMILES('C=[C]C(O)[CH]O[O]'),
    E0 = (300.778,'kJ/mol'),
    modes = [
        HarmonicOscillator(frequencies=([492.5,1135,1000,3615,1277.5,1000,3025,407.5,1350,352.5,2950,3100,1380,975,1025,1650,1380,1390,370,380,2900,435,1685,370,273.14,273.145],'cm^-1')),
        HinderedRotor(inertia=(0.0979926,'amu*angstrom^2'), symmetry=1, barrier=(5.1882,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.0979839,'amu*angstrom^2'), symmetry=1, barrier=(5.18828,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.0979968,'amu*angstrom^2'), symmetry=1, barrier=(5.18828,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.367914,'amu*angstrom^2'), symmetry=1, barrier=(19.4793,'kJ/mol'), semiclassical=False),
    ],
    spinMultiplicity = 4,
    opticalIsomers = 1,
    molecularWeight = (101.081,'amu'),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[0.869504,0.0742997,-0.00011163,9.25416e-08,-3.05441e-11,36282.8,29.4293], Tmin=(100,'K'), Tmax=(823.575,'K')), NASAPolynomial(coeffs=[9.32123,0.0265403,-1.24222e-05,2.34176e-09,-1.60341e-13,35118.2,-8.32324], Tmin=(823.575,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(300.778,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(266.063,'J/(mol*K)'), comment="""Thermo group additivity estimation: group(O2s-OsCs) + group(O2s-CsH) + group(O2s-OsH) + group(Cs-(Cds-Cds)CsOsH) + group(Cs-CsOsHH) + group(Cds-CdsCsH) + group(Cds-CdsHH) + radical(Cds_S) + radical(ROOJ) + radical(CCsJOOH)"""),
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
    E0 = (137.67,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS2',
    E0 = (245.53,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS3',
    E0 = (152.355,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS4',
    E0 = (290.733,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS5',
    E0 = (224.973,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS6',
    E0 = (137.67,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS7',
    E0 = (264.253,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS8',
    E0 = (137.67,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS9',
    E0 = (320.251,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS10',
    E0 = (295.706,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS11',
    E0 = (252.361,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS12',
    E0 = (326.083,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS13',
    E0 = (279.647,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS14',
    E0 = (279.647,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS15',
    E0 = (212.95,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS16',
    E0 = (247.967,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS17',
    E0 = (267.841,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS18',
    E0 = (192.101,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS19',
    E0 = (274.808,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS20',
    E0 = (384.282,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS21',
    E0 = (455.892,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS22',
    E0 = (586.93,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS23',
    E0 = (286.062,'kJ/mol'),
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
    E0 = (445.842,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS26',
    E0 = (554.721,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS27',
    E0 = (596.57,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS28',
    E0 = (215.917,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS29',
    E0 = (201.07,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS30',
    E0 = (226.638,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS31',
    E0 = (155.452,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS32',
    E0 = (137.67,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS33',
    E0 = (327.2,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS34',
    E0 = (220.568,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS35',
    E0 = (232.144,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS36',
    E0 = (277.039,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS37',
    E0 = (145.954,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS38',
    E0 = (145.954,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS39',
    E0 = (144.782,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS40',
    E0 = (380.303,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS41',
    E0 = (388.066,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS42',
    E0 = (647.45,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS43',
    E0 = (707.659,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

reaction(
    label = 'reaction1',
    reactants = ['C=[C]C(O)C([O])O[O](15278)'],
    products = ['[O]OC=O(5472)', 'C=C=CO(12571)'],
    transitionState = 'TS1',
    kinetics = Arrhenius(A=(5e+12,'s^-1'), n=0, Ea=(0,'kJ/mol'), T0=(1,'K'), Tmin=(300,'K'), Tmax=(1500,'K'), comment="""Exact match found for rate rule [RJJ]
Euclidian distance = 0
family: 1,4_Linear_birad_scission"""),
)

reaction(
    label = 'reaction2',
    reactants = ['H(8)', 'C=C=C(O)C([O])O[O](28905)'],
    products = ['C=[C]C(O)C([O])O[O](15278)'],
    transitionState = 'TS2',
    kinetics = Arrhenius(A=(169.619,'m^3/(mol*s)'), n=1.605, Ea=(12.4249,'kJ/mol'), T0=(1,'K'), comment="""Estimated using an average for rate rule [Cds_Ca;HJ]
Euclidian distance = 0
family: R_Addition_MultipleBond"""),
)

reaction(
    label = 'reaction3',
    reactants = ['H(8)', 'C=[C]C(O)C(=O)O[O](28906)'],
    products = ['C=[C]C(O)C([O])O[O](15278)'],
    transitionState = 'TS3',
    kinetics = Arrhenius(A=(1.83701,'m^3/(mol*s)'), n=1.71338, Ea=(46.5052,'kJ/mol'), T0=(1,'K'), comment="""Estimated using average of templates [CO_O;HJ] + [CO-NdNd_O;YJ] for rate rule [CO-NdNd_O;HJ]
Euclidian distance = 1.0
family: R_Addition_MultipleBond"""),
)

reaction(
    label = 'reaction4',
    reactants = ['H(8)', 'C#CC(O)C([O])O[O](28907)'],
    products = ['C=[C]C(O)C([O])O[O](15278)'],
    transitionState = 'TS4',
    kinetics = Arrhenius(A=(1.255e+11,'cm^3/(mol*s)'), n=1.005, Ea=(13.1503,'kJ/mol'), T0=(1,'K'), comment="""From training reaction 138 used for Ct-H_Ct-Cs;HJ
Exact match found for rate rule [Ct-H_Ct-Cs;HJ]
Euclidian distance = 0
family: R_Addition_MultipleBond"""),
)

reaction(
    label = 'reaction5',
    reactants = ['[O][CH]O[O](8201)', 'C=C=CO(12571)'],
    products = ['C=[C]C(O)C([O])O[O](15278)'],
    transitionState = 'TS5',
    kinetics = Arrhenius(A=(0.00401797,'m^3/(mol*s)'), n=2.41733, Ea=(22.1495,'kJ/mol'), T0=(1,'K'), comment="""Estimated using an average for rate rule [Cds_Ca;CJ]
Euclidian distance = 0
family: R_Addition_MultipleBond"""),
)

reaction(
    label = 'reaction6',
    reactants = ['[O]OC=O(5472)', '[CH2][C]=CO(18753)'],
    products = ['C=[C]C(O)C([O])O[O](15278)'],
    transitionState = 'TS6',
    kinetics = Arrhenius(A=(6.87291,'m^3/(mol*s)'), n=1.39198, Ea=(146.493,'kJ/mol'), T0=(1,'K'), comment="""Estimated using average of templates [CO_O;CJ] + [CO-NdH_O;YJ] for rate rule [CO-NdH_O;CJ]
Euclidian distance = 1.0
family: R_Addition_MultipleBond
Ea raised from 143.9 to 146.5 kJ/mol to match endothermicity of reaction."""),
)

reaction(
    label = 'reaction7',
    reactants = ['OH(D)(132)', 'C=C=CC([O])O[O](18038)'],
    products = ['C=[C]C(O)C([O])O[O](15278)'],
    transitionState = 'TS7',
    kinetics = Arrhenius(A=(986500,'cm^3/(mol*s)'), n=2.037, Ea=(0,'kJ/mol'), T0=(1,'K'), comment="""Estimated using template [Cds_Ca;OJ_pri] for rate rule [Cds-CsH_Ca;OJ_pri]
Euclidian distance = 1.0
family: R_Addition_MultipleBond
Ea raised from -6.0 to 0 kJ/mol."""),
)

reaction(
    label = 'reaction8',
    reactants = ['O2(2)', 'C=[C]C(O)C=O(27744)'],
    products = ['C=[C]C(O)C([O])O[O](15278)'],
    transitionState = 'TS8',
    kinetics = Arrhenius(A=(0.0849,'cm^3/(mol*s)'), n=3.486, Ea=(178.792,'kJ/mol'), T0=(1,'K'), Tmin=(300,'K'), Tmax=(1600,'K'), comment="""Estimated using template [CO-CsH_O;OJ] for rate rule [CO-CsH_O;O2b]
Euclidian distance = 1.0
Multiplied by reaction path degeneracy 2.0
family: R_Addition_MultipleBond
Ea raised from 176.4 to 178.8 kJ/mol to match endothermicity of reaction."""),
)

reaction(
    label = 'reaction9',
    reactants = ['C=[C]C(O)C([O])O[O](15278)'],
    products = ['[CH2]C=C(O)C([O])O[O](15274)'],
    transitionState = 'TS9',
    kinetics = Arrhenius(A=(3.677e+10,'s^-1'), n=0.839, Ea=(182.581,'kJ/mol'), T0=(1,'K'), comment="""Estimated using template [R2H_S;Cd_rad_out_Cd;Cs_H_out_noH] for rate rule [R2H_S;Cd_rad_out_Cd;Cs_H_out_NDMustO]
Euclidian distance = 2.0
family: intra_H_migration"""),
)

reaction(
    label = 'reaction10',
    reactants = ['C=[C]C(O)C([O])O[O](15278)'],
    products = ['C=[C]C(O)[C](O)O[O](28908)'],
    transitionState = 'TS10',
    kinetics = Arrhenius(A=(1.43381e+07,'s^-1'), n=1.70481, Ea=(158.036,'kJ/mol'), T0=(1,'K'), comment="""Estimated using template [R2H_S;Y_rad_out;Cs_H_out_NonDe] for rate rule [R2H_S;O_rad_out;Cs_H_out_NDMustO]
Euclidian distance = 1.41421356237
family: intra_H_migration"""),
)

reaction(
    label = 'reaction11',
    reactants = ['[CH]=CC(O)C([O])O[O](15282)'],
    products = ['C=[C]C(O)C([O])O[O](15278)'],
    transitionState = 'TS11',
    kinetics = Arrhenius(A=(1.08e+06,'s^-1'), n=1.99, Ea=(105.437,'kJ/mol'), T0=(1,'K'), comment="""From training reaction 17 used for R2H_D;Cd_rad_out_singleH;Cd_H_out_singleNd
Exact match found for rate rule [R2H_D;Cd_rad_out_singleH;Cd_H_out_singleNd]
Euclidian distance = 0
family: intra_H_migration"""),
)

reaction(
    label = 'reaction12',
    reactants = ['C=[C]C(O)[C]([O])OO(28909)'],
    products = ['C=[C]C(O)C([O])O[O](15278)'],
    transitionState = 'TS12',
    kinetics = Arrhenius(A=(40813.3,'s^-1'), n=2.17068, Ea=(135.172,'kJ/mol'), T0=(1,'K'), comment="""Estimated using average of templates [R3H_SS;Y_rad_out;O_H_out] + [R3H_SS_O;Y_rad_out;XH_out] for rate rule [R3H_SS_O;Y_rad_out;O_H_out]
Euclidian distance = 1.0
family: intra_H_migration"""),
)

reaction(
    label = 'reaction15',
    reactants = ['C=[C]C(O)C([O])O[O](15278)'],
    products = ['C=CC([O])C([O])O[O](12832)'],
    transitionState = 'TS13',
    kinetics = Arrhenius(A=(2.4115e+09,'s^-1'), n=1.00333, Ea=(141.977,'kJ/mol'), T0=(1,'K'), comment="""Estimated using template [R3H_SS_Cs;Cd_rad_out_Cd;XH_out] for rate rule [R3H_SS_Cs;Cd_rad_out_Cd;O_H_out]
Euclidian distance = 1.0
family: intra_H_migration"""),
)

reaction(
    label = 'reaction14',
    reactants = ['C=[C]C(O)C([O])O[O](15278)'],
    products = ['C=CC(O)[C]([O])O[O](15277)'],
    transitionState = 'TS14',
    kinetics = Arrhenius(A=(2.4115e+09,'s^-1'), n=1.00333, Ea=(141.977,'kJ/mol'), T0=(1,'K'), comment="""Estimated using an average for rate rule [R3H_SS_Cs;Cd_rad_out_Cd;XH_out]
Euclidian distance = 0
family: intra_H_migration"""),
)

reaction(
    label = 'reaction15',
    reactants = ['C=[C]C(O)C([O])O[O](15278)'],
    products = ['[CH2][C]=C(O)C(O)O[O](28910)'],
    transitionState = 'TS15',
    kinetics = Arrhenius(A=(111914,'s^-1'), n=2.27675, Ea=(75.2806,'kJ/mol'), T0=(1,'K'), comment="""Estimated using an average for rate rule [R3H_SS_Cs;O_rad_out;XH_out]
Euclidian distance = 0
family: intra_H_migration"""),
)

reaction(
    label = 'reaction16',
    reactants = ['C=[C]C([O])C(O)O[O](15280)'],
    products = ['C=[C]C(O)C([O])O[O](15278)'],
    transitionState = 'TS16',
    kinetics = Arrhenius(A=(94.0113,'s^-1'), n=2.81534, Ea=(105.641,'kJ/mol'), T0=(1,'K'), comment="""Estimated using average of templates [R4H_SSS;Y_rad_out;O_H_out] + [R4H_SSS;O_rad_out;XH_out] for rate rule [R4H_SSS;O_rad_out;O_H_out]
Euclidian distance = 1.0
family: intra_H_migration"""),
)

reaction(
    label = 'reaction17',
    reactants = ['C=[C]C(O)C([O])O[O](15278)'],
    products = ['[CH2][C]=C(O)C([O])OO(28911)'],
    transitionState = 'TS17',
    kinetics = Arrhenius(A=(1.17136e+07,'s^-1'), n=1.54267, Ea=(130.172,'kJ/mol'), T0=(1,'K'), comment="""Estimated using an average for rate rule [R4H_SSS_OCs;O_rad_out;XH_out]
Euclidian distance = 0
family: intra_H_migration"""),
)

reaction(
    label = 'reaction18',
    reactants = ['[CH]=[C]C(O)C(O)O[O](28912)'],
    products = ['C=[C]C(O)C([O])O[O](15278)'],
    transitionState = 'TS18',
    kinetics = Arrhenius(A=(136000,'s^-1'), n=1.9199, Ea=(33.0402,'kJ/mol'), T0=(1,'K'), comment="""Estimated using template [R5Hall;Cd_rad_out_singleH;XH_out] for rate rule [R5HJ_1;Cd_rad_out_singleH;O_H_out]
Euclidian distance = 1.41421356237
family: intra_H_migration"""),
)

reaction(
    label = 'reaction19',
    reactants = ['C=[C]C([O])C([O])OO(15284)'],
    products = ['C=[C]C(O)C([O])O[O](15278)'],
    transitionState = 'TS19',
    kinetics = Arrhenius(A=(1.37227e+06,'s^-1'), n=1.56745, Ea=(58.7826,'kJ/mol'), T0=(1,'K'), comment="""Estimated using template [R5H_SSSS;O_rad_out;XH_out] for rate rule [R5H_SSSS;O_rad_out;O_H_out]
Euclidian distance = 1.0
family: intra_H_migration"""),
)

reaction(
    label = 'reaction20',
    reactants = ['[CH]=[C]C(O)C([O])OO(28913)'],
    products = ['C=[C]C(O)C([O])O[O](15278)'],
    transitionState = 'TS20',
    kinetics = Arrhenius(A=(1.86943e+06,'s^-1'), n=1.85754, Ea=(151.521,'kJ/mol'), T0=(1,'K'), comment="""Estimated using template [RnH;Cd_rad_out_singleH;XH_out] for rate rule [R6HJ_1;Cd_rad_out_singleH;O_H_out]
Euclidian distance = 2.2360679775
family: intra_H_migration"""),
)

reaction(
    label = 'reaction21',
    reactants = ['OH(D)(132)', 'C=[C][CH]C([O])O[O](18045)'],
    products = ['C=[C]C(O)C([O])O[O](15278)'],
    transitionState = 'TS21',
    kinetics = Arrhenius(A=(3.05166e+07,'m^3/(mol*s)'), n=0.045, Ea=(0.1046,'kJ/mol'), T0=(1,'K'), comment="""Estimated using an average for rate rule [O_pri_rad;Y_rad]
Euclidian distance = 0
family: R_Recombination"""),
)

reaction(
    label = 'reaction22',
    reactants = ['H(8)', 'C=[C]C([O])C([O])O[O](15288)'],
    products = ['C=[C]C(O)C([O])O[O](15278)'],
    transitionState = 'TS22',
    kinetics = Arrhenius(A=(5.00518e+06,'m^3/(mol*s)'), n=0.282325, Ea=(7.09479,'kJ/mol'), T0=(1,'K'), comment="""Estimated using average of templates [Y_rad;O_rad/NonDe] + [H_rad;O_sec_rad] for rate rule [H_rad;O_rad/NonDe]
Euclidian distance = 1.0
family: R_Recombination"""),
)

reaction(
    label = 'reaction23',
    reactants = ['O2(2)', 'C=[C]C(O)[CH][O](28102)'],
    products = ['C=[C]C(O)C([O])O[O](15278)'],
    transitionState = 'TS23',
    kinetics = Arrhenius(A=(2.18266e+06,'m^3/(mol*s)'), n=0.193158, Ea=(0,'kJ/mol'), T0=(1,'K'), comment="""Estimated using an average for rate rule [Y_rad;O2_birad]
Euclidian distance = 0
Multiplied by reaction path degeneracy 2.0
family: R_Recombination
Ea raised from -25.1 to 0 kJ/mol."""),
)

reaction(
    label = 'reaction24',
    reactants = ['[O][CH]O[O](8201)', '[CH2][C]=CO(18753)'],
    products = ['C=[C]C(O)C([O])O[O](15278)'],
    transitionState = 'TS24',
    kinetics = Arrhenius(A=(1.9789e+07,'m^3/(mol*s)'), n=-0.126319, Ea=(0,'kJ/mol'), T0=(1,'K'), comment="""Estimated using an average for rate rule [Y_rad;Y_rad]
Euclidian distance = 0
family: R_Recombination
Ea raised from -15.6 to -15.6 kJ/mol.
Ea raised from -15.6 to 0 kJ/mol."""),
)

reaction(
    label = 'reaction25',
    reactants = ['H(8)', '[CH2][C]=C(O)C([O])O[O](28914)'],
    products = ['C=[C]C(O)C([O])O[O](15278)'],
    transitionState = 'TS25',
    kinetics = Arrhenius(A=(4.34078e+06,'m^3/(mol*s)'), n=0.278577, Ea=(0,'kJ/mol'), T0=(1,'K'), comment="""Estimated using an average for rate rule [Y_rad;H_rad]
Euclidian distance = 0
family: R_Recombination
Ea raised from -1.4 to 0 kJ/mol."""),
)

reaction(
    label = 'reaction26',
    reactants = ['H(8)', 'C=[C]C(O)[C]([O])O[O](28915)'],
    products = ['C=[C]C(O)C([O])O[O](15278)'],
    transitionState = 'TS26',
    kinetics = Arrhenius(A=(4.34078e+06,'m^3/(mol*s)'), n=0.278577, Ea=(0,'kJ/mol'), T0=(1,'K'), comment="""Estimated using an average for rate rule [Y_rad;H_rad]
Euclidian distance = 0
family: R_Recombination
Ea raised from -1.4 to 0 kJ/mol."""),
)

reaction(
    label = 'reaction27',
    reactants = ['H(8)', '[CH]=[C]C(O)C([O])O[O](28916)'],
    products = ['C=[C]C(O)C([O])O[O](15278)'],
    transitionState = 'TS27',
    kinetics = Arrhenius(A=(4.34078e+06,'m^3/(mol*s)'), n=0.278577, Ea=(0,'kJ/mol'), T0=(1,'K'), comment="""Estimated using an average for rate rule [Y_rad;H_rad]
Euclidian distance = 0
family: R_Recombination
Ea raised from -1.4 to 0 kJ/mol."""),
)

reaction(
    label = 'reaction28',
    reactants = ['C=[C]C(O)C([O])O[O](15278)'],
    products = ['C=C=C(O)C(O)O[O](28917)'],
    transitionState = 'TS28',
    kinetics = Arrhenius(A=(2.00399e+09,'s^-1'), n=0.37, Ea=(78.2471,'kJ/mol'), T0=(1,'K'), comment="""Estimated using average of templates [R3;Y_rad;XH_Rrad_De] + [R3radExo;Y_rad;XH_Rrad] for rate rule [R3radExo;Y_rad;XH_Rrad_De]
Euclidian distance = 1.0
family: Intra_Disproportionation"""),
)

reaction(
    label = 'reaction29',
    reactants = ['C=[C]C(O)C([O])O[O](15278)'],
    products = ['C=[C]C(O)C(=O)OO(27760)'],
    transitionState = 'TS29',
    kinetics = Arrhenius(A=(7.437e+08,'s^-1'), n=1.045, Ea=(63.4002,'kJ/mol'), T0=(1,'K'), comment="""Estimated using an average for rate rule [R3radExo;Y_rad;XH_Rrad]
Euclidian distance = 0
family: Intra_Disproportionation"""),
)

reaction(
    label = 'reaction30',
    reactants = ['C=[C]C(O)C([O])O[O](15278)'],
    products = ['C=CC(O)C(=O)O[O](15293)'],
    transitionState = 'TS30',
    kinetics = Arrhenius(A=(2.6374e+09,'s^-1'), n=0.37, Ea=(88.9686,'kJ/mol'), T0=(1,'K'), comment="""Estimated using average of templates [R3;Y_rad_De;XH_Rrad] + [R3radExo;Y_rad;XH_Rrad] for rate rule [R3radExo;Y_rad_De;XH_Rrad]
Euclidian distance = 1.0
family: Intra_Disproportionation"""),
)

reaction(
    label = 'reaction31',
    reactants = ['C=[C]C(O)C([O])O[O](15278)'],
    products = ['C=C=C(O)C([O])OO(28918)'],
    transitionState = 'TS31',
    kinetics = Arrhenius(A=(2.59e+08,'s^-1'), n=0.311, Ea=(17.782,'kJ/mol'), T0=(1,'K'), Tmin=(300,'K'), Tmax=(1500,'K'), comment="""Estimated using template [R4;Y_rad;XH_Rrad_De] for rate rule [R4radExo;Y_rad;XH_Rrad_De]
Euclidian distance = 1.0
family: Intra_Disproportionation"""),
)

reaction(
    label = 'reaction32',
    reactants = ['C=[C]C(O)C([O])O[O](15278)'],
    products = ['O2(S)(5486)', 'C=[C]C(O)C=O(27744)'],
    transitionState = 'TS32',
    kinetics = Arrhenius(A=(5e+12,'s^-1'), n=0, Ea=(0,'kJ/mol'), T0=(1,'K'), Tmin=(300,'K'), Tmax=(1500,'K'), comment="""Exact match found for rate rule [RJJ]
Euclidian distance = 0
family: 1,4_Linear_birad_scission"""),
)

reaction(
    label = 'reaction33',
    reactants = ['C=[C]C(O)C([O])O[O](15278)'],
    products = ['O(T)(63)', 'C=[C]C(O)C1OO1(28919)'],
    transitionState = 'TS33',
    kinetics = Arrhenius(A=(2.30161e+09,'s^-1'), n=1.08, Ea=(189.53,'kJ/mol'), T0=(1,'K'), comment="""Estimated using an average for rate rule [R2OO_S;Y_rad_intra;OOJ]
Euclidian distance = 0
family: Cyclic_Ether_Formation"""),
)

reaction(
    label = 'reaction34',
    reactants = ['C=[C]C(O)C([O])O[O](15278)'],
    products = ['O(T)(63)', 'C=C1OC([O])C1O(28920)'],
    transitionState = 'TS34',
    kinetics = Arrhenius(A=(3.11355e+11,'s^-1'), n=0, Ea=(82.8981,'kJ/mol'), T0=(1,'K'), comment="""Estimated using template [R3OO_SS;Y_rad_intra;OO] for rate rule [R3OO_SS;Cd_rad_out;OOJ]
Euclidian distance = 1.41421356237
family: Cyclic_Ether_Formation"""),
)

reaction(
    label = 'reaction35',
    reactants = ['C=[C]C(O)C([O])O[O](15278)'],
    products = ['[CH2]C(=CO)C([O])O[O](15193)'],
    transitionState = 'TS35',
    kinetics = Arrhenius(A=(8.66e+11,'s^-1'), n=0.438, Ea=(94.4747,'kJ/mol'), T0=(1,'K'), comment="""Estimated using an average for rate rule [cCs(-HR!H)CJ;CdsJ;C]
Euclidian distance = 0
family: 1,2_shiftC"""),
)

reaction(
    label = 'reaction36',
    reactants = ['C=[C]C(O)C([O])O[O](15278)'],
    products = ['[O]O(16)', 'C=[C]C(O)=C[O](28921)'],
    transitionState = 'TS36',
    kinetics = Arrhenius(A=(9.58174e+10,'s^-1'), n=0.573333, Ea=(139.369,'kJ/mol'), T0=(1,'K'), comment="""Estimated using an average for rate rule [R2OO]
Euclidian distance = 0
family: HO2_Elimination_from_PeroxyRadical"""),
)

reaction(
    label = 'reaction37',
    reactants = ['C=[C]C(O)C([O])O[O](15278)'],
    products = ['C=[C]C(O)C1OOO1(28922)'],
    transitionState = 'TS37',
    kinetics = Arrhenius(A=(1.62e+12,'s^-1'), n=-0.305, Ea=(8.28432,'kJ/mol'), T0=(1,'K'), comment="""Estimated using template [R4_SSS;Y_rad_out;Ypri_rad_out] for rate rule [R4_SSS;O_rad;Opri_rad]
Euclidian distance = 1.41421356237
family: Birad_recombination"""),
)

reaction(
    label = 'reaction38',
    reactants = ['C=[C]C(O)C([O])O[O](15278)'],
    products = ['C=C1OC(O[O])C1O(28884)'],
    transitionState = 'TS38',
    kinetics = Arrhenius(A=(1.62e+12,'s^-1'), n=-0.305, Ea=(8.28432,'kJ/mol'), T0=(1,'K'), comment="""Estimated using template [R4_SSS;Y_rad_out;Ypri_rad_out] for rate rule [R4_SSS;Y_rad_out;Opri_rad]
Euclidian distance = 1.0
family: Birad_recombination"""),
)

reaction(
    label = 'reaction39',
    reactants = ['C=[C]C(O)C([O])O[O](15278)'],
    products = ['C=C1OOC([O])C1O(28923)'],
    transitionState = 'TS39',
    kinetics = Arrhenius(A=(7.76e+09,'s^-1'), n=0.311, Ea=(7.1128,'kJ/mol'), T0=(1,'K'), comment="""Estimated using template [R5_SSSS;Y_rad_out;Ypri_rad_out] for rate rule [R5_SSSS;O_rad;Ypri_rad_out]
Euclidian distance = 1.0
family: Birad_recombination"""),
)

reaction(
    label = 'reaction40',
    reactants = ['C=[C][CH]C(O[O])OO(28924)'],
    products = ['C=[C]C(O)C([O])O[O](15278)'],
    transitionState = 'TS40',
    kinetics = Arrhenius(A=(3.01978e+11,'s^-1'), n=0, Ea=(107.111,'kJ/mol'), T0=(1,'K'), comment="""Estimated using an average for rate rule [R2OOH_S;Y_rad_out]
Euclidian distance = 0
family: intra_OH_migration"""),
)

reaction(
    label = 'reaction41',
    reactants = ['O(T)(63)', 'C=[C]C(O)C([O])[O](28925)'],
    products = ['C=[C]C(O)C([O])O[O](15278)'],
    transitionState = 'TS41',
    kinetics = Arrhenius(A=(87.1677,'m^3/(mol*s)'), n=1.88017, Ea=(5.1666,'kJ/mol'), T0=(1,'K'), Tmin=(303.03,'K'), Tmax=(2000,'K'), comment="""From training reaction 3 used for O_rad/NonDe;O_birad
Exact match found for rate rule [O_rad/NonDe;O_birad]
Euclidian distance = 0
Multiplied by reaction path degeneracy 2.0
family: Birad_R_Recombination"""),
)

reaction(
    label = 'reaction42',
    reactants = ['[C]=C(584)', '[O]OC([O])[CH]O(8411)'],
    products = ['C=[C]C(O)C([O])O[O](15278)'],
    transitionState = 'TS42',
    kinetics = Arrhenius(A=(1.14854e+06,'m^3/(mol*s)'), n=0.575199, Ea=(34.3157,'kJ/mol'), T0=(1,'K'), comment="""Estimated using template [Y_rad;Birad] for rate rule [C_rad/H/CsO;Birad]
Euclidian distance = 4.0
family: Birad_R_Recombination"""),
)

reaction(
    label = 'reaction43',
    reactants = ['O(T)(63)', 'C=[C]C(O)[CH]O[O](28926)'],
    products = ['C=[C]C(O)C([O])O[O](15278)'],
    transitionState = 'TS43',
    kinetics = Arrhenius(A=(93609.6,'m^3/(mol*s)'), n=1.13083, Ea=(163.847,'kJ/mol'), T0=(1,'K'), Tmin=(303.03,'K'), Tmax=(2000,'K'), comment="""Estimated using template [Y_rad;O_birad] for rate rule [C_rad/H/CsO;O_birad]
Euclidian distance = 4.0
family: Birad_R_Recombination"""),
)

network(
    label = '5159',
    isomers = [
        'C=[C]C(O)C([O])O[O](15278)',
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
    label = '5159',
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

