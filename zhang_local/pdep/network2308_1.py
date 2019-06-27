species(
    label = '[O]O[CH]OC([O])O[O](8194)',
    structure = SMILES('[O]O[CH]OC([O])O[O]'),
    E0 = (91.4629,'kJ/mol'),
    modes = [
        HarmonicOscillator(frequencies=([470,515,1100,1170,900,1100,1380,1390,370,380,2900,435,3025,407.5,1350,352.5,180,180,1878.3,1878.71],'cm^-1')),
        HinderedRotor(inertia=(0.217455,'amu*angstrom^2'), symmetry=1, barrier=(4.99973,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.217814,'amu*angstrom^2'), symmetry=1, barrier=(5.00798,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.217315,'amu*angstrom^2'), symmetry=1, barrier=(4.99649,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.217098,'amu*angstrom^2'), symmetry=1, barrier=(4.99151,'kJ/mol'), semiclassical=False),
    ],
    spinMultiplicity = 5,
    opticalIsomers = 1,
    molecularWeight = (122.034,'amu'),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[0.764865,0.0980434,-0.000230046,2.49431e-07,-9.5886e-11,11091.2,29.8539], Tmin=(100,'K'), Tmax=(870.413,'K')), NASAPolynomial(coeffs=[-2.6804,0.0469213,-2.65619e-05,5.2032e-09,-3.55719e-13,14227.3,60.5668], Tmin=(870.413,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(91.4629,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(216.176,'J/(mol*K)'), comment="""Thermo group additivity estimation: group(O2s-CsCs) + group(O2s-OsCs) + group(O2s-OsCs) + group(O2s-CsH) + group(O2s-OsH) + group(O2s-OsH) + group(Cs-OsOsOsH) + group(Cs-OsOsHH) + radical(OCOJ) + radical(OCJO) + radical(ROOJ) + radical(ROOJ)"""),
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
    label = '[O]O[CH]OC(=O)O[O](8886)',
    structure = SMILES('[O]O[CH]OC(=O)O[O]'),
    E0 = (-194.712,'kJ/mol'),
    modes = [
        HarmonicOscillator(frequencies=([470,515,1100,1170,900,1100,3025,407.5,1350,352.5,200,800,960,1120,1280,1440,1600],'cm^-1')),
        HinderedRotor(inertia=(0.156089,'amu*angstrom^2'), symmetry=1, barrier=(3.5888,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.156089,'amu*angstrom^2'), symmetry=1, barrier=(3.5888,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.156089,'amu*angstrom^2'), symmetry=1, barrier=(3.5888,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.156089,'amu*angstrom^2'), symmetry=1, barrier=(3.5888,'kJ/mol'), semiclassical=False),
    ],
    spinMultiplicity = 4,
    opticalIsomers = 1,
    molecularWeight = (121.026,'amu'),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[0.24497,0.070285,-8.39888e-05,4.52416e-08,-9.0631e-12,-23272.8,30.3107], Tmin=(100,'K'), Tmax=(1346.8,'K')), NASAPolynomial(coeffs=[21.6935,0.000330376,8.8695e-07,-2.19003e-10,1.53783e-14,-28483.1,-77.4466], Tmin=(1346.8,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(-194.712,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(191.233,'J/(mol*K)'), comment="""Thermo group additivity estimation: group(O2s-Cs(Cds-O2d)) + group(O2s-OsCs) + group(O2s-O2s(Cds-O2d)) + group(O2s-OsH) + group(O2s-OsH) + group(Cs-OsOsHH) + group(Cds-OdOsOs) + radical(C(=O)OOJ) + radical(ROOJ) + radical(OCJO)"""),
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
    label = '[O]OC([O])OC=O(8887)',
    structure = SMILES('[O]OC([O])OC=O'),
    E0 = (-300.843,'kJ/mol'),
    modes = [
        HarmonicOscillator(frequencies=([492.5,1135,1000,1380,1390,370,380,2900,435,2782.5,750,1395,475,1775,1000,180,1234.59,1234.8],'cm^-1')),
        HinderedRotor(inertia=(0.13179,'amu*angstrom^2'), symmetry=1, barrier=(3.03011,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.00279597,'amu*angstrom^2'), symmetry=1, barrier=(3.02122,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.131415,'amu*angstrom^2'), symmetry=1, barrier=(3.02148,'kJ/mol'), semiclassical=False),
    ],
    spinMultiplicity = 3,
    opticalIsomers = 1,
    molecularWeight = (106.034,'amu'),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[2.31933,0.0542553,-0.00012052,1.46402e-07,-6.26597e-11,-36139.4,23.4189], Tmin=(100,'K'), Tmax=(830.597,'K')), NASAPolynomial(coeffs=[-5.8208,0.0479697,-2.70228e-05,5.42372e-09,-3.81646e-13,-33218.2,70.6248], Tmin=(830.597,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(-300.843,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(195.39,'J/(mol*K)'), comment="""Thermo group additivity estimation: group(O2s-Cs(Cds-O2d)) + group(O2s-OsCs) + group(O2s-CsH) + group(O2s-OsH) + group(Cs-OsOsOsH) + group(Cds-OdOsH) + radical(OCOJ) + radical(ROOJ)"""),
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
    label = '[O]O[CH]OC=O(8157)',
    structure = SMILES('[O]O[CH]OC=O'),
    E0 = (-149.557,'kJ/mol'),
    modes = [
        HarmonicOscillator(frequencies=([492.5,1135,1000,2782.5,750,1395,475,1775,1000,3025,407.5,1350,352.5,296.346,296.397],'cm^-1')),
        HinderedRotor(inertia=(0.00325132,'amu*angstrom^2'), symmetry=1, barrier=(18.8438,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.593368,'amu*angstrom^2'), symmetry=1, barrier=(36.965,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.944592,'amu*angstrom^2'), symmetry=1, barrier=(58.9036,'kJ/mol'), semiclassical=False),
    ],
    spinMultiplicity = 3,
    opticalIsomers = 1,
    molecularWeight = (90.0349,'amu'),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[2.11148,0.0439904,-5.22323e-05,3.28477e-08,-8.38967e-12,-17921.6,19.7044], Tmin=(100,'K'), Tmax=(945.088,'K')), NASAPolynomial(coeffs=[8.71534,0.0160391,-7.86758e-06,1.55157e-09,-1.10726e-13,-19169.8,-11.7823], Tmin=(945.088,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(-149.557,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(170.447,'J/(mol*K)'), comment="""Thermo group additivity estimation: group(O2s-Cs(Cds-O2d)) + group(O2s-OsCs) + group(O2s-OsH) + group(Cs-OsOsHH) + group(Cds-OdOsH) + radical(ROOJ) + radical(OCJO)"""),
)

species(
    label = '[O]O[CH]O[C](O)O[O](8888)',
    structure = SMILES('[O]O[CH]O[C](O)O[O]'),
    E0 = (69.3114,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
    molecularWeight = (122.034,'amu'),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[0.240084,0.101688,-0.000216159,2.1228e-07,-7.59337e-11,8453.58,31.9533], Tmin=(100,'K'), Tmax=(882.644,'K')), NASAPolynomial(coeffs=[6.18068,0.0294269,-1.63035e-05,3.12968e-09,-2.09927e-13,9171.01,14.0396], Tmin=(882.644,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(69.3114,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(212.019,'J/(mol*K)'), comment="""Thermo group additivity estimation: group(O2s-CsCs) + group(O2s-OsCs) + group(O2s-OsCs) + group(O2s-CsH) + group(O2s-OsH) + group(O2s-OsH) + group(Cs-OsOsOsH) + group(Cs-OsOsHH) + radical(ROOJ) + radical(OCJO) + radical(Cs_P) + radical(ROOJ)"""),
)

species(
    label = '[O]OCO[C]([O])O[O](8889)',
    structure = SMILES('[O]OCO[C]([O])O[O]'),
    E0 = (107.749,'kJ/mol'),
    modes = [
        HarmonicOscillator(frequencies=([470,515,1100,1170,900,1100,360,370,350,2750,2850,1437.5,1250,1305,750,350,180,180,1186,1190.63],'cm^-1')),
        HinderedRotor(inertia=(0.209887,'amu*angstrom^2'), symmetry=1, barrier=(4.82571,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.209938,'amu*angstrom^2'), symmetry=1, barrier=(4.82689,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.211976,'amu*angstrom^2'), symmetry=1, barrier=(4.87374,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.210279,'amu*angstrom^2'), symmetry=1, barrier=(4.83472,'kJ/mol'), semiclassical=False),
    ],
    spinMultiplicity = 5,
    opticalIsomers = 1,
    molecularWeight = (122.034,'amu'),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[0.864471,0.0822566,-0.000159317,1.57024e-07,-5.81654e-11,13059.3,30.3953], Tmin=(100,'K'), Tmax=(841.717,'K')), NASAPolynomial(coeffs=[5.16051,0.0324063,-1.80256e-05,3.5712e-09,-2.48421e-13,13378.8,16.6036], Tmin=(841.717,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(107.749,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(216.176,'J/(mol*K)'), comment="""Thermo group additivity estimation: group(O2s-CsCs) + group(O2s-OsCs) + group(O2s-OsCs) + group(O2s-CsH) + group(O2s-OsH) + group(O2s-OsH) + group(Cs-OsOsOsH) + group(Cs-OsOsHH) + radical(ROOJ) + radical(Cs_P) + radical(OCOJ) + radical(ROOJ)"""),
)

species(
    label = '[O]O[CH]O[C]([O])OO(8890)',
    structure = SMILES('[O]O[CH]O[C]([O])OO'),
    E0 = (144.705,'kJ/mol'),
    modes = [
        HarmonicOscillator(frequencies=([3025,407.5,1350,352.5,3615,1310,387.5,850,1000,360,370,350,492.5,1135,1000,200,800,1200,1600],'cm^-1')),
        HinderedRotor(inertia=(0.156089,'amu*angstrom^2'), symmetry=1, barrier=(3.5888,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.156089,'amu*angstrom^2'), symmetry=1, barrier=(3.5888,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.156089,'amu*angstrom^2'), symmetry=1, barrier=(3.5888,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.156089,'amu*angstrom^2'), symmetry=1, barrier=(3.5888,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.156089,'amu*angstrom^2'), symmetry=1, barrier=(3.5888,'kJ/mol'), semiclassical=False),
    ],
    spinMultiplicity = 5,
    opticalIsomers = 1,
    molecularWeight = (122.034,'amu'),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[0.529444,0.102485,-0.000236283,2.53193e-07,-9.70283e-11,17503.8,30.6271], Tmin=(100,'K'), Tmax=(864.156,'K')), NASAPolynomial(coeffs=[-1.20359,0.0461209,-2.66862e-05,5.27766e-09,-3.63249e-13,20207.4,52.6448], Tmin=(864.156,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(144.705,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(212.019,'J/(mol*K)'), comment="""Thermo group additivity estimation: group(O2s-CsCs) + group(O2s-OsCs) + group(O2s-OsCs) + group(O2s-CsH) + group(O2s-OsH) + group(O2s-OsH) + group(Cs-OsOsOsH) + group(Cs-OsOsHH) + radical(OCOJ) + radical(Cs_P) + radical(ROOJ) + radical(OCJO)"""),
)

species(
    label = '[O]O[C]([O])O[CH]OO(8891)',
    structure = SMILES('[O]O[C]([O])O[CH]OO'),
    E0 = (144.705,'kJ/mol'),
    modes = [
        HarmonicOscillator(frequencies=([3025,407.5,1350,352.5,3615,1310,387.5,850,1000,360,370,350,492.5,1135,1000,200,800,1200,1600],'cm^-1')),
        HinderedRotor(inertia=(0.156089,'amu*angstrom^2'), symmetry=1, barrier=(3.5888,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.156089,'amu*angstrom^2'), symmetry=1, barrier=(3.5888,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.156089,'amu*angstrom^2'), symmetry=1, barrier=(3.5888,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.156089,'amu*angstrom^2'), symmetry=1, barrier=(3.5888,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.156089,'amu*angstrom^2'), symmetry=1, barrier=(3.5888,'kJ/mol'), semiclassical=False),
    ],
    spinMultiplicity = 5,
    opticalIsomers = 1,
    molecularWeight = (122.034,'amu'),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[0.529444,0.102485,-0.000236283,2.53193e-07,-9.70283e-11,17503.8,30.6271], Tmin=(100,'K'), Tmax=(864.156,'K')), NASAPolynomial(coeffs=[-1.20359,0.0461209,-2.66862e-05,5.27766e-09,-3.63249e-13,20207.4,52.6448], Tmin=(864.156,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(144.705,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(212.019,'J/(mol*K)'), comment="""Thermo group additivity estimation: group(O2s-CsCs) + group(O2s-OsCs) + group(O2s-OsCs) + group(O2s-CsH) + group(O2s-OsH) + group(O2s-OsH) + group(Cs-OsOsOsH) + group(Cs-OsOsHH) + radical(OCJO) + radical(OCOJ) + radical(ROOJ) + radical(Cs_P)"""),
)

species(
    label = '[O][CH]O[CH]O[O](8404)',
    structure = SMILES('[O][CH]O[CH]O[O]'),
    E0 = (246.637,'kJ/mol'),
    modes = [
        HarmonicOscillator(frequencies=([492.5,1135,1000,3000,3050,390,425,1340,1360,335,370,180,180,987.086,987.583],'cm^-1')),
        HinderedRotor(inertia=(0.233233,'amu*angstrom^2'), symmetry=1, barrier=(5.36249,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.233045,'amu*angstrom^2'), symmetry=1, barrier=(5.35816,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.2333,'amu*angstrom^2'), symmetry=1, barrier=(5.36402,'kJ/mol'), semiclassical=False),
    ],
    spinMultiplicity = 5,
    opticalIsomers = 1,
    molecularWeight = (90.0349,'amu'),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[1.34999,0.0738578,-0.000159407,1.6125e-07,-5.90184e-11,29744.2,22.2321], Tmin=(100,'K'), Tmax=(880.549,'K')), NASAPolynomial(coeffs=[3.43888,0.026574,-1.44775e-05,2.77816e-09,-1.87011e-13,30841.6,20.7402], Tmin=(880.549,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(246.637,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(170.447,'J/(mol*K)'), comment="""Thermo group additivity estimation: group(O2s-CsCs) + group(O2s-OsCs) + group(O2s-CsH) + group(O2s-OsH) + group(Cs-OsOsHH) + group(Cs-OsOsHH) + radical(OCOJ) + radical(OCJO) + radical(OCJO) + radical(ROOJ)"""),
)

species(
    label = '[O]O[CH]O[C]([O])O[O](8892)',
    structure = SMILES('[O]O[CH]O[C]([O])O[O]'),
    E0 = (296.709,'kJ/mol'),
    modes = [
        HarmonicOscillator(frequencies=([470,515,1100,1170,900,1100,360,370,350,3025,407.5,1350,352.5,180,1728.98,1729.07,1729.9],'cm^-1')),
        HinderedRotor(inertia=(0.275994,'amu*angstrom^2'), symmetry=1, barrier=(6.34565,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.275706,'amu*angstrom^2'), symmetry=1, barrier=(6.33903,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.275783,'amu*angstrom^2'), symmetry=1, barrier=(6.34078,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.276074,'amu*angstrom^2'), symmetry=1, barrier=(6.34749,'kJ/mol'), semiclassical=False),
    ],
    spinMultiplicity = 6,
    opticalIsomers = 1,
    molecularWeight = (121.026,'amu'),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[0.799655,0.0988009,-0.000240549,2.6168e-07,-1.00277e-10,35774,30.4172], Tmin=(100,'K'), Tmax=(874.263,'K')), NASAPolynomial(coeffs=[-2.85975,0.0445799,-2.57658e-05,5.05396e-09,-3.44553e-13,39125.9,63.0902], Tmin=(874.263,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(296.709,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(191.233,'J/(mol*K)'), comment="""Thermo group additivity estimation: group(O2s-CsCs) + group(O2s-OsCs) + group(O2s-OsCs) + group(O2s-CsH) + group(O2s-OsH) + group(O2s-OsH) + group(Cs-OsOsOsH) + group(Cs-OsOsHH) + radical(OCJO) + radical(OCOJ) + radical(ROOJ) + radical(Cs_P) + radical(ROOJ)"""),
)

species(
    label = '[O]OCOC(=O)O[O](8893)',
    structure = SMILES('[O]OCOC(=O)O[O]'),
    E0 = (-383.673,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
    molecularWeight = (122.034,'amu'),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[0.666881,0.0495747,1.16135e-05,-7.78072e-08,4.08206e-11,-46003,29.0047], Tmin=(100,'K'), Tmax=(919.357,'K')), NASAPolynomial(coeffs=[28.6294,-0.0100221,7.58787e-06,-1.45831e-09,9.14547e-14,-53767.4,-117.812], Tmin=(919.357,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(-383.673,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(216.176,'J/(mol*K)'), comment="""Thermo group additivity estimation: group(O2s-Cs(Cds-O2d)) + group(O2s-OsCs) + group(O2s-O2s(Cds-O2d)) + group(O2s-OsH) + group(O2s-OsH) + group(Cs-OsOsHH) + group(Cds-OdOsOs) + radical(C(=O)OOJ) + radical(ROOJ)"""),
)

species(
    label = '[O]O[CH]OC(=O)OO(8894)',
    structure = SMILES('[O]O[CH]OC(=O)OO'),
    E0 = (-389.101,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
    molecularWeight = (122.034,'amu'),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[-0.201905,0.0759093,-8.58154e-05,4.35945e-08,-8.22692e-12,-46632.7,31.1638], Tmin=(100,'K'), Tmax=(1444.19,'K')), NASAPolynomial(coeffs=[24.0983,0.000619013,6.78527e-07,-1.615e-10,1.03369e-14,-52818.7,-92.1191], Tmin=(1444.19,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(-389.101,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(212.019,'J/(mol*K)'), comment="""Thermo group additivity estimation: group(O2s-Cs(Cds-O2d)) + group(O2s-OsCs) + group(O2s-O2s(Cds-O2d)) + group(O2s-OsH) + group(O2s-OsH) + group(Cs-OsOsHH) + group(Cds-OdOsOs) + radical(OCJO) + radical(ROOJ)"""),
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
    label = '[O]O[CH]OC1OO1(8895)',
    structure = SMILES('[O]O[CH]OC1OO1'),
    E0 = (34.7663,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
    molecularWeight = (106.034,'amu'),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[1.29508,0.0620123,-9.24744e-05,7.06095e-08,-2.10212e-11,4276.44,23.5103], Tmin=(100,'K'), Tmax=(866.819,'K')), NASAPolynomial(coeffs=[11.0901,0.0144682,-6.14405e-06,1.09311e-09,-7.19981e-14,2666.42,-21.837], Tmin=(866.819,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(34.7663,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(195.39,'J/(mol*K)'), comment="""Thermo group additivity estimation: group(O2s-CsCs) + group(O2s-OsCs) + group(O2s-OsCs) + group(O2s-OsCs) + group(O2s-OsH) + group(Cs-OsOsOsH) + group(Cs-OsOsHH) + ring(dioxirane) + radical(ROOJ) + radical(OCJO)"""),
)

species(
    label = '[O]OC1OC([O])O1(8896)',
    structure = SMILES('[O]OC1OC([O])O1'),
    E0 = (-161.968,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
    molecularWeight = (106.034,'amu'),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[5.79444,0.0107871,1.29696e-05,-1.24797e-08,2.42711e-12,-19591,9.02057], Tmin=(100,'K'), Tmax=(2128.33,'K')), NASAPolynomial(coeffs=[23.1158,0.0110849,-1.03934e-05,2.09069e-09,-1.36252e-13,-34404.8,-105.108], Tmin=(2128.33,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(-161.968,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(203.705,'J/(mol*K)'), comment="""Thermo group additivity estimation: group(O2s-CsCs) + group(O2s-CsCs) + group(O2s-OsCs) + group(O2s-CsH) + group(O2s-OsH) + group(Cs-OsOsOsH) + group(Cs-OsOsOsH) + ring(Cyclobutane) + radical(ROOJ) + radical(OCOJ)"""),
)

species(
    label = '[O]OC([O])OC1OO1(8897)',
    structure = SMILES('[O]OC([O])OC1OO1'),
    E0 = (-116.52,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
    molecularWeight = (122.034,'amu'),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[1.69681,0.0698051,-0.000151136,1.70001e-07,-6.83551e-11,-13949.7,26.5408], Tmin=(100,'K'), Tmax=(851.525,'K')), NASAPolynomial(coeffs=[-3.7116,0.0468697,-2.55788e-05,5.03268e-09,-3.48595e-13,-11276,62.0549], Tmin=(851.525,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(-116.52,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(220.334,'J/(mol*K)'), comment="""Thermo group additivity estimation: group(O2s-CsCs) + group(O2s-OsCs) + group(O2s-OsCs) + group(O2s-OsCs) + group(O2s-CsH) + group(O2s-OsH) + group(Cs-OsOsOsH) + group(Cs-OsOsOsH) + ring(dioxirane) + radical(OCOJ) + radical(ROOJ)"""),
)

species(
    label = '[O]OC1OC(O[O])O1(8197)',
    structure = SMILES('[O]OC1OC(O[O])O1'),
    E0 = (-165.856,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
    molecularWeight = (122.034,'amu'),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[4.71907,0.0297531,-1.36069e-05,1.32607e-09,9.1894e-14,-20020.3,12.4158], Tmin=(100,'K'), Tmax=(2694.31,'K')), NASAPolynomial(coeffs=[47.83,-0.0168504,2.6519e-06,-3.00076e-10,2.03815e-14,-49566.5,-250.018], Tmin=(2694.31,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(-165.856,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(224.491,'J/(mol*K)'), comment="""Thermo group additivity estimation: group(O2s-CsCs) + group(O2s-CsCs) + group(O2s-OsCs) + group(O2s-OsCs) + group(O2s-OsH) + group(O2s-OsH) + group(Cs-OsOsOsH) + group(Cs-OsOsOsH) + ring(Cyclobutane) + radical(ROOJ) + radical(ROOJ)"""),
)

species(
    label = '[O]O[CH]OC1OOO1(8898)',
    structure = SMILES('[O]O[CH]OC1OOO1'),
    E0 = (75.9788,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
    molecularWeight = (122.034,'amu'),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[1.53959,0.0540272,-4.96842e-05,2.18842e-08,-3.86349e-12,9226.54,25.7278], Tmin=(100,'K'), Tmax=(1336.6,'K')), NASAPolynomial(coeffs=[13.406,0.0185147,-9.8301e-06,2.00582e-09,-1.45386e-13,6054.41,-34.9637], Tmin=(1336.6,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(75.9788,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(220.334,'J/(mol*K)'), comment="""Thermo group additivity estimation: group(O2s-CsCs) + group(O2s-OsCs) + group(O2s-OsCs) + group(O2s-OsCs) + group(O2s-OsOs) + group(O2s-OsH) + group(Cs-OsOsOsH) + group(Cs-OsOsHH) + ring(Cyclobutane) + radical(ROOJ) + radical(OCJO)"""),
)

species(
    label = '[O]OC1OOC([O])O1(8899)',
    structure = SMILES('[O]OC1OOC([O])O1'),
    E0 = (-179.669,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
    molecularWeight = (122.034,'amu'),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[3.86206,0.028874,-1.15436e-05,1.46518e-09,-1.5973e-14,-21631.6,17.2012], Tmin=(100,'K'), Tmax=(2825.16,'K')), NASAPolynomial(coeffs=[34.3568,-0.00418596,6.38419e-07,-1.42041e-10,1.40963e-14,-42899.2,-168.734], Tmin=(2825.16,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(-179.669,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(228.648,'J/(mol*K)'), comment="""Thermo group additivity estimation: group(O2s-CsCs) + group(O2s-OsCs) + group(O2s-OsCs) + group(O2s-OsCs) + group(O2s-CsH) + group(O2s-OsH) + group(Cs-OsOsOsH) + group(Cs-OsOsOsH) + ring(124trioxolane) + radical(OCOJ) + radical(ROOJ)"""),
)

species(
    label = '[CH]O[O](2819)',
    structure = SMILES('[CH]O[O]'),
    E0 = (465.475,'kJ/mol'),
    modes = [
        HarmonicOscillator(frequencies=([492.5,1135,1000,744.599,4000],'cm^-1')),
        HinderedRotor(inertia=(0.274125,'amu*angstrom^2'), symmetry=1, barrier=(6.30268,'kJ/mol'), semiclassical=False),
    ],
    spinMultiplicity = 4,
    opticalIsomers = 1,
    molecularWeight = (45.0174,'amu'),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[3.27108,0.0184967,-3.44448e-05,3.15738e-08,-1.07949e-11,56007.6,9.95572], Tmin=(100,'K'), Tmax=(890.352,'K')), NASAPolynomial(coeffs=[4.92998,0.00531001,-2.56884e-06,4.73001e-10,-3.11601e-14,55939.5,3.42146], Tmin=(890.352,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(465.475,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(78.9875,'J/(mol*K)'), comment="""Thermo group additivity estimation: group(O2s-OsCs) + group(O2s-OsH) + group(Cs-OsHHH) + radical(ROOJ) + radical(CH2_triplet)"""),
)

species(
    label = '[O]OC([O])[O](4172)',
    structure = SMILES('[O]OC([O])[O]'),
    E0 = (77.6018,'kJ/mol'),
    modes = [
        HarmonicOscillator(frequencies=([492.5,1135,1000,1380,1390,370,380,2900,435,2220.26,2220.69],'cm^-1')),
        HinderedRotor(inertia=(3.41843e-05,'amu*angstrom^2'), symmetry=1, barrier=(0.119627,'kJ/mol'), semiclassical=False),
    ],
    spinMultiplicity = 4,
    opticalIsomers = 1,
    molecularWeight = (77.0162,'amu'),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[3.3674,0.0423788,-0.000145172,1.9796e-07,-8.61253e-11,9328.67,15.7326], Tmin=(100,'K'), Tmax=(861.03,'K')), NASAPolynomial(coeffs=[-15.5083,0.0534406,-3.09507e-05,6.16376e-09,-4.27232e-13,15419.6,120.467], Tmin=(861.03,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(77.6018,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(128.874,'J/(mol*K)'), comment="""Thermo group additivity estimation: group(O2s-OsCs) + group(O2s-CsH) + group(O2s-CsH) + group(O2s-OsH) + group(Cs-OsOsOsH) + radical(OCOJ) + radical(ROOJ) + radical(OCOJ)"""),
)

species(
    label = '[O]O[CH]OC([O])[O](8900)',
    structure = SMILES('[O]O[CH]OC([O])[O]'),
    E0 = (95.3511,'kJ/mol'),
    modes = [
        HarmonicOscillator(frequencies=([492.5,1135,1000,1380,1390,370,380,2900,435,3025,407.5,1350,352.5,180,1717.13,1719.27,1719.85,1721.89],'cm^-1')),
        HinderedRotor(inertia=(0.183853,'amu*angstrom^2'), symmetry=1, barrier=(4.22714,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.184607,'amu*angstrom^2'), symmetry=1, barrier=(4.24448,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.183352,'amu*angstrom^2'), symmetry=1, barrier=(4.21562,'kJ/mol'), semiclassical=False),
    ],
    spinMultiplicity = 5,
    opticalIsomers = 1,
    molecularWeight = (106.034,'amu'),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[1.75086,0.0816534,-0.000218032,2.60497e-07,-1.06231e-10,11518.1,24.573], Tmin=(100,'K'), Tmax=(863.793,'K')), NASAPolynomial(coeffs=[-11.3321,0.0589225,-3.38814e-05,6.71037e-09,-4.62994e-13,16886.5,103.766], Tmin=(863.793,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(95.3511,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(195.39,'J/(mol*K)'), comment="""Thermo group additivity estimation: group(O2s-CsCs) + group(O2s-OsCs) + group(O2s-CsH) + group(O2s-CsH) + group(O2s-OsH) + group(Cs-OsOsOsH) + group(Cs-OsOsHH) + radical(OCJO) + radical(OCOJ) + radical(OCOJ) + radical(ROOJ)"""),
)

species(
    label = '[O][CH]OC([O])O[O](8901)',
    structure = SMILES('[O][CH]OC([O])O[O]'),
    E0 = (95.3511,'kJ/mol'),
    modes = [
        HarmonicOscillator(frequencies=([492.5,1135,1000,1380,1390,370,380,2900,435,3025,407.5,1350,352.5,180,1718.06,1718.28,1720.26,1721.31],'cm^-1')),
        HinderedRotor(inertia=(0.18422,'amu*angstrom^2'), symmetry=1, barrier=(4.23557,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.183596,'amu*angstrom^2'), symmetry=1, barrier=(4.22123,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.183909,'amu*angstrom^2'), symmetry=1, barrier=(4.22842,'kJ/mol'), semiclassical=False),
    ],
    spinMultiplicity = 5,
    opticalIsomers = 1,
    molecularWeight = (106.034,'amu'),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[1.75065,0.0816561,-0.000218043,2.60512e-07,-1.06239e-10,11518.2,25.2668], Tmin=(100,'K'), Tmax=(863.78,'K')), NASAPolynomial(coeffs=[-11.3316,0.0589217,-3.38809e-05,6.71025e-09,-4.62984e-13,16886.3,104.457], Tmin=(863.78,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(95.3511,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(195.39,'J/(mol*K)'), comment="""Thermo group additivity estimation: group(O2s-CsCs) + group(O2s-OsCs) + group(O2s-CsH) + group(O2s-CsH) + group(O2s-OsH) + group(Cs-OsOsOsH) + group(Cs-OsOsHH) + radical(OCOJ) + radical(ROOJ) + radical(OCJO) + radical(OCOJ)"""),
)

species(
    label = '[CH]OC([O])O[O](8400)',
    structure = SMILES('[CH]OC([O])O[O]'),
    E0 = (331.938,'kJ/mol'),
    modes = [
        HarmonicOscillator(frequencies=([492.5,1135,1000,1380,1390,370,380,2900,435,180,180,2106.73,2106.98,2108.16,2108.18],'cm^-1')),
        HinderedRotor(inertia=(0.129957,'amu*angstrom^2'), symmetry=1, barrier=(2.98798,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.131221,'amu*angstrom^2'), symmetry=1, barrier=(3.01702,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.129936,'amu*angstrom^2'), symmetry=1, barrier=(2.98749,'kJ/mol'), semiclassical=False),
    ],
    spinMultiplicity = 5,
    opticalIsomers = 1,
    molecularWeight = (90.0349,'amu'),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[2.05365,0.0655934,-0.000166056,1.93586e-07,-7.82481e-11,39971.1,21.8361], Tmin=(100,'K'), Tmax=(859.228,'K')), NASAPolynomial(coeffs=[-5.67417,0.0431564,-2.49126e-05,4.95398e-09,-3.43084e-13,43455.3,70.4934], Tmin=(859.228,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(331.938,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(170.447,'J/(mol*K)'), comment="""Thermo group additivity estimation: group(O2s-CsCs) + group(O2s-OsCs) + group(O2s-CsH) + group(O2s-OsH) + group(Cs-OsOsOsH) + group(Cs-OsHHH) + radical(ROOJ) + radical(OCOJ) + radical(CH2_triplet)"""),
)

species(
    label = '[O]O[CH]O[CH]O[O](8902)',
    structure = SMILES('[O]O[CH]O[CH]O[O]'),
    E0 = (242.749,'kJ/mol'),
    modes = [
        HarmonicOscillator(frequencies=([470,515,1100,1170,900,1100,3000,3050,390,425,1340,1360,335,370,258.116,258.148,258.226],'cm^-1')),
        HinderedRotor(inertia=(0.171711,'amu*angstrom^2'), symmetry=1, barrier=(8.12311,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.171677,'amu*angstrom^2'), symmetry=1, barrier=(8.12303,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.171684,'amu*angstrom^2'), symmetry=1, barrier=(8.12269,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.336557,'amu*angstrom^2'), symmetry=1, barrier=(15.9268,'kJ/mol'), semiclassical=False),
    ],
    spinMultiplicity = 5,
    opticalIsomers = 1,
    molecularWeight = (106.034,'amu'),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[0.346899,0.0904809,-0.000172411,1.51779e-07,-4.95267e-11,29318,26.1862], Tmin=(100,'K'), Tmax=(909.545,'K')), NASAPolynomial(coeffs=[12.072,0.0146016,-7.17361e-06,1.27448e-09,-8.00138e-14,28190.8,-23.7402], Tmin=(909.545,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(242.749,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(191.233,'J/(mol*K)'), comment="""Thermo group additivity estimation: group(O2s-CsCs) + group(O2s-OsCs) + group(O2s-OsCs) + group(O2s-OsH) + group(O2s-OsH) + group(Cs-OsOsHH) + group(Cs-OsOsHH) + radical(OCJO) + radical(ROOJ) + radical(ROOJ) + radical(OCJO)"""),
)

species(
    label = '[O]O[C]OC([O])O[O](8903)',
    structure = SMILES('[O]O[C]OC([O])O[O]'),
    E0 = (363.491,'kJ/mol'),
    modes = [
        HarmonicOscillator(frequencies=([470,515,1100,1170,900,1100,1380,1390,370,380,2900,435,200,800,1066.67,1333.33,1600],'cm^-1')),
        HinderedRotor(inertia=(0.156089,'amu*angstrom^2'), symmetry=1, barrier=(3.5888,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.156089,'amu*angstrom^2'), symmetry=1, barrier=(3.5888,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.156089,'amu*angstrom^2'), symmetry=1, barrier=(3.5888,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.156089,'amu*angstrom^2'), symmetry=1, barrier=(3.5888,'kJ/mol'), semiclassical=False),
    ],
    spinMultiplicity = 6,
    opticalIsomers = 1,
    molecularWeight = (121.026,'amu'),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[0.619928,0.0873951,-0.000173689,1.67451e-07,-6.06341e-11,43827,29.3175], Tmin=(100,'K'), Tmax=(841.638,'K')), NASAPolynomial(coeffs=[7.64349,0.0263412,-1.5556e-05,3.12512e-09,-2.18209e-13,43624.9,2.46649], Tmin=(841.638,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(363.491,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(191.233,'J/(mol*K)'), comment="""Thermo group additivity estimation: group(O2s-CsCs) + group(O2s-OsCs) + group(O2s-OsCs) + group(O2s-CsH) + group(O2s-OsH) + group(O2s-OsH) + group(Cs-OsOsOsH) + group(Cs-OsOsHH) + radical(OCOJ) + radical(ROOJ) + radical(ROOJ) + radical(CH2_triplet)"""),
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
    E0 = (91.4629,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS2',
    E0 = (91.4629,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS3',
    E0 = (91.5007,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS4',
    E0 = (91.4629,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS5',
    E0 = (91.4629,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS6',
    E0 = (249.499,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS7',
    E0 = (261.929,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS8',
    E0 = (279.877,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS9',
    E0 = (196.836,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS10',
    E0 = (457.776,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS11',
    E0 = (238.01,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS12',
    E0 = (508.514,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS13',
    E0 = (154.863,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS14',
    E0 = (154.863,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS15',
    E0 = (91.4629,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS16',
    E0 = (279.493,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS17',
    E0 = (197.764,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS18',
    E0 = (99.3233,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS19',
    E0 = (99.7472,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS20',
    E0 = (99.7472,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS21',
    E0 = (98.5757,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS22',
    E0 = (548.243,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS23',
    E0 = (343.552,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS24',
    E0 = (502.232,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS25',
    E0 = (357.627,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS26',
    E0 = (649.63,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS27',
    E0 = (575.296,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

reaction(
    label = 'reaction1',
    reactants = ['[O]O[CH]OC([O])O[O](8194)'],
    products = ['[O]OC=O(5472)', '[O]OC=O(5472)'],
    transitionState = 'TS1',
    kinetics = Arrhenius(A=(5e+12,'s^-1'), n=0, Ea=(0,'kJ/mol'), T0=(1,'K'), Tmin=(300,'K'), Tmax=(1500,'K'), comment="""Exact match found for rate rule [RJJ]
Euclidian distance = 0
family: 1,4_Linear_birad_scission"""),
)

reaction(
    label = 'reaction2',
    reactants = ['H(8)', '[O]O[CH]OC(=O)O[O](8886)'],
    products = ['[O]O[CH]OC([O])O[O](8194)'],
    transitionState = 'TS2',
    kinetics = Arrhenius(A=(1.83701,'m^3/(mol*s)'), n=1.71338, Ea=(74.3707,'kJ/mol'), T0=(1,'K'), comment="""Estimated using average of templates [CO_O;HJ] + [CO-NdNd_O;YJ] for rate rule [CO-NdNd_O;HJ]
Euclidian distance = 1.0
family: R_Addition_MultipleBond
Ea raised from 71.6 to 74.4 kJ/mol to match endothermicity of reaction."""),
)

reaction(
    label = 'reaction3',
    reactants = ['O(T)(63)', '[O]OC([O])OC=O(8887)'],
    products = ['[O]O[CH]OC([O])O[O](8194)'],
    transitionState = 'TS3',
    kinetics = Arrhenius(A=(4000,'m^3/(mol*s)'), n=1.39, Ea=(149.31,'kJ/mol'), T0=(1,'K'), comment="""Estimated using template [Od_CO-NdH;YJ] for rate rule [Od_CO-NdH;O_atom_triplet]
Euclidian distance = 2.0
family: R_Addition_MultipleBond"""),
)

reaction(
    label = 'reaction4',
    reactants = ['[O]OC=O(5472)', '[O][CH]O[O](8201)'],
    products = ['[O]O[CH]OC([O])O[O](8194)'],
    transitionState = 'TS4',
    kinetics = Arrhenius(A=(4000,'m^3/(mol*s)'), n=1.39, Ea=(58.0705,'kJ/mol'), T0=(1,'K'), comment="""Estimated using template [Od_CO-NdH;YJ] for rate rule [Od_CO-NdH;CJ]
Euclidian distance = 1.0
family: R_Addition_MultipleBond
Ea raised from 56.5 to 58.1 kJ/mol to match endothermicity of reaction."""),
)

reaction(
    label = 'reaction5',
    reactants = ['O2(2)', '[O]O[CH]OC=O(8157)'],
    products = ['[O]O[CH]OC([O])O[O](8194)'],
    transitionState = 'TS5',
    kinetics = Arrhenius(A=(8.49e-08,'m^3/(mol*s)'), n=3.486, Ea=(249.647,'kJ/mol'), T0=(1,'K'), comment="""Estimated using template [CO-NdH_O;OJ] for rate rule [CO-NdH_O;O2b]
Euclidian distance = 1.0
Multiplied by reaction path degeneracy 2.0
family: R_Addition_MultipleBond
Ea raised from 248.0 to 249.6 kJ/mol to match endothermicity of reaction."""),
)

reaction(
    label = 'reaction6',
    reactants = ['[O]O[CH]OC([O])O[O](8194)'],
    products = ['[O]O[CH]O[C](O)O[O](8888)'],
    transitionState = 'TS6',
    kinetics = Arrhenius(A=(1.43381e+07,'s^-1'), n=1.70481, Ea=(158.036,'kJ/mol'), T0=(1,'K'), comment="""Estimated using template [R2H_S;Y_rad_out;Cs_H_out_NonDe] for rate rule [R2H_S;O_rad_out;Cs_H_out_NDMustO]
Euclidian distance = 1.41421356237
family: intra_H_migration"""),
)

reaction(
    label = 'reaction7',
    reactants = ['[O]OCO[C]([O])O[O](8889)'],
    products = ['[O]O[CH]OC([O])O[O](8194)'],
    transitionState = 'TS7',
    kinetics = Arrhenius(A=(6e+08,'s^-1'), n=1.23, Ea=(154.18,'kJ/mol'), T0=(1,'K'), Tmin=(300,'K'), Tmax=(1500,'K'), comment="""Estimated using an average for rate rule [R3H_SS_O;Y_rad_out;Cs_H_out_H/NonDeO]
Euclidian distance = 0
Multiplied by reaction path degeneracy 2.0
family: intra_H_migration"""),
)

reaction(
    label = 'reaction8',
    reactants = ['[O]O[CH]O[C]([O])OO(8890)'],
    products = ['[O]O[CH]OC([O])O[O](8194)'],
    transitionState = 'TS8',
    kinetics = Arrhenius(A=(40813.3,'s^-1'), n=2.17068, Ea=(135.172,'kJ/mol'), T0=(1,'K'), comment="""Estimated using average of templates [R3H_SS;Y_rad_out;O_H_out] + [R3H_SS_O;Y_rad_out;XH_out] for rate rule [R3H_SS_O;Y_rad_out;O_H_out]
Euclidian distance = 1.0
family: intra_H_migration"""),
)

reaction(
    label = 'reaction9',
    reactants = ['[O]O[C]([O])O[CH]OO(8891)'],
    products = ['[O]O[CH]OC([O])O[O](8194)'],
    transitionState = 'TS9',
    kinetics = Arrhenius(A=(1.10205e+06,'s^-1'), n=1.54368, Ea=(52.1315,'kJ/mol'), T0=(1,'K'), comment="""Estimated using template [R5Hall;Y_rad_out;XH_out] for rate rule [R5HJ_2;Y_rad_out;O_H_out]
Euclidian distance = 1.41421356237
family: intra_H_migration"""),
)

reaction(
    label = 'reaction10',
    reactants = ['[O][CH]O[O](8201)', '[O][CH]O[O](8201)'],
    products = ['[O]O[CH]OC([O])O[O](8194)'],
    transitionState = 'TS10',
    kinetics = Arrhenius(A=(1.9789e+07,'m^3/(mol*s)'), n=-0.126319, Ea=(0,'kJ/mol'), T0=(1,'K'), comment="""Estimated using an average for rate rule [Y_rad;Y_rad]
Euclidian distance = 0
family: R_Recombination
Ea raised from -15.6 to -15.6 kJ/mol.
Ea raised from -15.6 to 0 kJ/mol."""),
)

reaction(
    label = 'reaction11',
    reactants = ['O2(2)', '[O][CH]O[CH]O[O](8404)'],
    products = ['[O]O[CH]OC([O])O[O](8194)'],
    transitionState = 'TS11',
    kinetics = Arrhenius(A=(2.18266e+06,'m^3/(mol*s)'), n=0.193158, Ea=(0,'kJ/mol'), T0=(1,'K'), comment="""Estimated using an average for rate rule [Y_rad;O2_birad]
Euclidian distance = 0
Multiplied by reaction path degeneracy 2.0
family: R_Recombination
Ea raised from -25.1 to 0 kJ/mol."""),
)

reaction(
    label = 'reaction12',
    reactants = ['H(8)', '[O]O[CH]O[C]([O])O[O](8892)'],
    products = ['[O]O[CH]OC([O])O[O](8194)'],
    transitionState = 'TS12',
    kinetics = Arrhenius(A=(4.34078e+06,'m^3/(mol*s)'), n=0.278577, Ea=(0,'kJ/mol'), T0=(1,'K'), comment="""Estimated using an average for rate rule [Y_rad;H_rad]
Euclidian distance = 0
family: R_Recombination
Ea raised from -1.4 to 0 kJ/mol."""),
)

reaction(
    label = 'reaction13',
    reactants = ['[O]O[CH]OC([O])O[O](8194)'],
    products = ['[O]OCOC(=O)O[O](8893)'],
    transitionState = 'TS13',
    kinetics = Arrhenius(A=(7.437e+08,'s^-1'), n=1.045, Ea=(63.4002,'kJ/mol'), T0=(1,'K'), comment="""Estimated using an average for rate rule [R3radExo;Y_rad_NDe;XH_Rrad]
Euclidian distance = 0
family: Intra_Disproportionation"""),
)

reaction(
    label = 'reaction14',
    reactants = ['[O]O[CH]OC([O])O[O](8194)'],
    products = ['[O]O[CH]OC(=O)OO(8894)'],
    transitionState = 'TS14',
    kinetics = Arrhenius(A=(7.437e+08,'s^-1'), n=1.045, Ea=(63.4002,'kJ/mol'), T0=(1,'K'), comment="""Estimated using an average for rate rule [R3radExo;Y_rad;XH_Rrad]
Euclidian distance = 0
family: Intra_Disproportionation"""),
)

reaction(
    label = 'reaction15',
    reactants = ['[O]O[CH]OC([O])O[O](8194)'],
    products = ['O2(S)(5486)', '[O]O[CH]OC=O(8157)'],
    transitionState = 'TS15',
    kinetics = Arrhenius(A=(5e+12,'s^-1'), n=0, Ea=(0,'kJ/mol'), T0=(1,'K'), Tmin=(300,'K'), Tmax=(1500,'K'), comment="""Exact match found for rate rule [RJJ]
Euclidian distance = 0
family: 1,4_Linear_birad_scission"""),
)

reaction(
    label = 'reaction16',
    reactants = ['[O]O[CH]OC([O])O[O](8194)'],
    products = ['O(T)(63)', '[O]O[CH]OC1OO1(8895)'],
    transitionState = 'TS16',
    kinetics = Arrhenius(A=(2.30161e+09,'s^-1'), n=1.08, Ea=(188.03,'kJ/mol'), T0=(1,'K'), comment="""Estimated using an average for rate rule [R2OO_S;Y_rad_intra;OOJ]
Euclidian distance = 0
family: Cyclic_Ether_Formation"""),
)

reaction(
    label = 'reaction17',
    reactants = ['[O]O[CH]OC([O])O[O](8194)'],
    products = ['O(T)(63)', '[O]OC1OC([O])O1(8896)'],
    transitionState = 'TS17',
    kinetics = Arrhenius(A=(2.58279e+10,'s^-1'), n=0.53, Ea=(106.301,'kJ/mol'), T0=(1,'K'), comment="""Estimated using average of templates [RnOO;C_sec_rad_intra;OOJ] + [R3OO;C_sec_rad_intra;OO] for rate rule [R3OO;C_rad/H/NonDeO_intra;OOJ]
Euclidian distance = 1.41421356237
family: Cyclic_Ether_Formation"""),
)

reaction(
    label = 'reaction18',
    reactants = ['[O]O[CH]OC([O])O[O](8194)'],
    products = ['[O]OC([O])OC1OO1(8897)'],
    transitionState = 'TS18',
    kinetics = Arrhenius(A=(1.61967e+11,'s^-1'), n=0.0247333, Ea=(7.86034,'kJ/mol'), T0=(1,'K'), comment="""Estimated using template [Rn;Y_rad_out;Cpri_rad_out_single] for rate rule [R3_SS;O_rad;Cpri_rad_out_H/NonDeO]
Euclidian distance = 3.0
family: Birad_recombination"""),
)

reaction(
    label = 'reaction19',
    reactants = ['[O]O[CH]OC([O])O[O](8194)'],
    products = ['[O]OC1OC(O[O])O1(8197)'],
    transitionState = 'TS19',
    kinetics = Arrhenius(A=(1.62e+12,'s^-1'), n=-0.305, Ea=(8.28432,'kJ/mol'), T0=(1,'K'), comment="""Estimated using template [R4_SSS;C_rad_out_single;Ypri_rad_out] for rate rule [R4_SSS;C_rad_out_H/NonDeO;Opri_rad]
Euclidian distance = 2.2360679775
family: Birad_recombination"""),
)

reaction(
    label = 'reaction20',
    reactants = ['[O]O[CH]OC([O])O[O](8194)'],
    products = ['[O]O[CH]OC1OOO1(8898)'],
    transitionState = 'TS20',
    kinetics = Arrhenius(A=(1.62e+12,'s^-1'), n=-0.305, Ea=(8.28432,'kJ/mol'), T0=(1,'K'), comment="""Estimated using template [R4_SSS;Y_rad_out;Ypri_rad_out] for rate rule [R4_SSS;O_rad;Opri_rad]
Euclidian distance = 1.41421356237
family: Birad_recombination"""),
)

reaction(
    label = 'reaction21',
    reactants = ['[O]O[CH]OC([O])O[O](8194)'],
    products = ['[O]OC1OOC([O])O1(8899)'],
    transitionState = 'TS21',
    kinetics = Arrhenius(A=(7.76e+09,'s^-1'), n=0.311, Ea=(7.1128,'kJ/mol'), T0=(1,'K'), comment="""Estimated using template [R5_SSSS;C_rad_out_single;Ypri_rad_out] for rate rule [R5_SSSS;C_rad_out_H/NonDeO;Opri_rad]
Euclidian distance = 2.2360679775
family: Birad_recombination"""),
)

reaction(
    label = 'reaction22',
    reactants = ['[CH]O[O](2819)', '[O]OC([O])[O](4172)'],
    products = ['[O]O[CH]OC([O])O[O](8194)'],
    transitionState = 'TS22',
    kinetics = Arrhenius(A=(87.1677,'m^3/(mol*s)'), n=1.88017, Ea=(5.1666,'kJ/mol'), T0=(1,'K'), Tmin=(303.03,'K'), Tmax=(2000,'K'), comment="""Estimated using an average for rate rule [O_rad/NonDe;Birad]
Euclidian distance = 0
Multiplied by reaction path degeneracy 2.0
family: Birad_R_Recombination"""),
)

reaction(
    label = 'reaction23',
    reactants = ['O(T)(63)', '[O]O[CH]OC([O])[O](8900)'],
    products = ['[O]O[CH]OC([O])O[O](8194)'],
    transitionState = 'TS23',
    kinetics = Arrhenius(A=(87.1677,'m^3/(mol*s)'), n=1.88017, Ea=(5.1666,'kJ/mol'), T0=(1,'K'), Tmin=(303.03,'K'), Tmax=(2000,'K'), comment="""From training reaction 3 used for O_rad/NonDe;O_birad
Exact match found for rate rule [O_rad/NonDe;O_birad]
Euclidian distance = 0
Multiplied by reaction path degeneracy 2.0
family: Birad_R_Recombination"""),
)

reaction(
    label = 'reaction24',
    reactants = ['O(T)(63)', '[O][CH]OC([O])O[O](8901)'],
    products = ['[O]O[CH]OC([O])O[O](8194)'],
    transitionState = 'TS24',
    kinetics = Arrhenius(A=(93609.6,'m^3/(mol*s)'), n=1.13083, Ea=(163.847,'kJ/mol'), T0=(1,'K'), Tmin=(303.03,'K'), Tmax=(2000,'K'), comment="""From training reaction 2 used for Y_rad;O_birad
Exact match found for rate rule [Y_rad;O_birad]
Euclidian distance = 0
family: Birad_R_Recombination"""),
)

reaction(
    label = 'reaction25',
    reactants = ['O2(2)', '[CH]OC([O])O[O](8400)'],
    products = ['[O]O[CH]OC([O])O[O](8194)'],
    transitionState = 'TS25',
    kinetics = Arrhenius(A=(2.29708e+06,'m^3/(mol*s)'), n=0.575199, Ea=(34.3157,'kJ/mol'), T0=(1,'K'), comment="""Estimated using an average for rate rule [Y_rad;Birad]
Euclidian distance = 0
Multiplied by reaction path degeneracy 2.0
family: Birad_R_Recombination"""),
)

reaction(
    label = 'reaction26',
    reactants = ['O(T)(63)', '[O]O[CH]O[CH]O[O](8902)'],
    products = ['[O]O[CH]OC([O])O[O](8194)'],
    transitionState = 'TS26',
    kinetics = Arrhenius(A=(187219,'m^3/(mol*s)'), n=1.13083, Ea=(163.847,'kJ/mol'), T0=(1,'K'), Tmin=(303.03,'K'), Tmax=(2000,'K'), comment="""Estimated using template [Y_rad;O_birad] for rate rule [C_rad/H/O2;O_birad]
Euclidian distance = 4.0
Multiplied by reaction path degeneracy 2.0
family: Birad_R_Recombination"""),
)

reaction(
    label = 'reaction27',
    reactants = ['H(8)', '[O]O[C]OC([O])O[O](8903)'],
    products = ['[O]O[CH]OC([O])O[O](8194)'],
    transitionState = 'TS27',
    kinetics = Arrhenius(A=(1e+07,'m^3/(mol*s)'), n=0, Ea=(0,'kJ/mol'), T0=(1,'K'), comment="""Estimated using an average for rate rule [H_rad;Birad]
Euclidian distance = 0
family: Birad_R_Recombination"""),
)

network(
    label = '2308',
    isomers = [
        '[O]O[CH]OC([O])O[O](8194)',
    ],
    reactants = [
        ('[O]OC=O(5472)', '[O]OC=O(5472)'),
    ],
    bathGas = {
        'N2': 0.25,
        'Ne': 0.25,
        'He': 0.25,
        'Ar': 0.25,
    },
)

pressureDependence(
    label = '2308',
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

