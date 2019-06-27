species(
    label = '[O]OC([O])C([O])O[O](8193)',
    structure = SMILES('[O]OC([O])C([O])O[O]'),
    E0 = (81.5691,'kJ/mol'),
    modes = [
        HarmonicOscillator(frequencies=([470,515,1100,1170,900,1100,1380,1383.33,1386.67,1390,370,373.333,376.667,380,2800,3000,430,440,180,180,2538.39],'cm^-1')),
        HinderedRotor(inertia=(0.154837,'amu*angstrom^2'), symmetry=1, barrier=(3.56001,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.154625,'amu*angstrom^2'), symmetry=1, barrier=(3.55514,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.154343,'amu*angstrom^2'), symmetry=1, barrier=(3.54864,'kJ/mol'), semiclassical=False),
    ],
    spinMultiplicity = 5,
    opticalIsomers = 1,
    molecularWeight = (122.034,'amu'),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[0.918447,0.086913,-0.00018522,1.92506e-07,-7.26533e-11,9903.11,31.9668], Tmin=(100,'K'), Tmax=(864.789,'K')), NASAPolynomial(coeffs=[1.30575,0.0387657,-2.13022e-05,4.15652e-09,-2.84744e-13,11569.5,40.1766], Tmin=(864.789,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(81.5691,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(220.334,'J/(mol*K)'), comment="""Thermo group additivity estimation: group(O2s-OsCs) + group(O2s-OsCs) + group(O2s-CsH) + group(O2s-CsH) + group(O2s-OsH) + group(O2s-OsH) + group(Cs-CsOsOsH) + group(Cs-CsOsOsH) + radical(ROOJ) + radical(ROOJ) + radical(CCOJ) + radical(CCOJ)"""),
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
    label = '[O]OC(=O)C([O])O[O](8904)',
    structure = SMILES('[O]OC(=O)C([O])O[O]'),
    E0 = (-144.027,'kJ/mol'),
    modes = [
        HarmonicOscillator(frequencies=([470,515,1100,1170,900,1100,1380,1390,370,380,2900,435,180,252.432,643.245,643.257,643.274,2657.4],'cm^-1')),
        HinderedRotor(inertia=(0.24815,'amu*angstrom^2'), symmetry=1, barrier=(72.8609,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.0948389,'amu*angstrom^2'), symmetry=1, barrier=(2.18053,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(3.16897,'amu*angstrom^2'), symmetry=1, barrier=(72.8608,'kJ/mol'), semiclassical=False),
    ],
    spinMultiplicity = 4,
    opticalIsomers = 1,
    molecularWeight = (121.026,'amu'),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[1.39871,0.0582439,-7.14177e-05,4.4105e-08,-1.07465e-11,-17229.6,29.8291], Tmin=(100,'K'), Tmax=(1003.04,'K')), NASAPolynomial(coeffs=[12.1776,0.015259,-7.13561e-06,1.37999e-09,-9.76292e-14,-19391.9,-22.2053], Tmin=(1003.04,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(-144.027,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(195.39,'J/(mol*K)'), comment="""Thermo group additivity estimation: group(O2s-OsCs) + group(O2s-O2s(Cds-O2d)) + group(O2s-CsH) + group(O2s-OsH) + group(O2s-OsH) + group(Cs-CsOsOsH) + group(Cds-OdCsOs) + radical(C=OCOJ) + radical(ROOJ) + radical(C(=O)OOJ)"""),
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
    label = '[O]OC([O])C=O(8158)',
    structure = SMILES('[O]OC([O])C=O'),
    E0 = (-70.5683,'kJ/mol'),
    modes = [
        HarmonicOscillator(frequencies=([492.5,1135,1000,1380,1390,370,380,2900,435,2782.5,750,1395,475,1775,1000,2561.73],'cm^-1')),
        HinderedRotor(inertia=(0.384019,'amu*angstrom^2'), symmetry=1, barrier=(8.82935,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.383579,'amu*angstrom^2'), symmetry=1, barrier=(8.81924,'kJ/mol'), semiclassical=False),
    ],
    spinMultiplicity = 3,
    opticalIsomers = 1,
    molecularWeight = (90.0349,'amu'),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[2.15222,0.0453291,-7.21521e-05,6.47296e-08,-2.29111e-11,-8425.39,22.8205], Tmin=(100,'K'), Tmax=(824.463,'K')), NASAPolynomial(coeffs=[5.8206,0.0192119,-9.49892e-06,1.82843e-09,-1.2643e-13,-8747.53,7.54555], Tmin=(824.463,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(-70.5683,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(174.604,'J/(mol*K)'), comment="""Thermo group additivity estimation: group(O2s-OsCs) + group(O2s-CsH) + group(O2s-OsH) + group(Cs-CsOsOsH) + group(Cds-OdCsH) + radical(C=OCOJ) + radical(ROOJ)"""),
)

species(
    label = '[O]O[C](O)C([O])O[O](8905)',
    structure = SMILES('[O]O[C](O)C([O])O[O]'),
    E0 = (61.1103,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
    molecularWeight = (122.034,'amu'),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[0.607375,0.0902088,-0.000182617,1.78466e-07,-6.43439e-11,7457.15,33.6051], Tmin=(100,'K'), Tmax=(870.231,'K')), NASAPolynomial(coeffs=[5.80398,0.0300313,-1.63353e-05,3.15924e-09,-2.14721e-13,7926.87,17.1523], Tmin=(870.231,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(61.1103,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(216.176,'J/(mol*K)'), comment="""Thermo group additivity estimation: group(O2s-OsCs) + group(O2s-OsCs) + group(O2s-CsH) + group(O2s-CsH) + group(O2s-OsH) + group(O2s-OsH) + group(Cs-CsOsOsH) + group(Cs-CsOsOsH) + radical(Cs_P) + radical(CCOJ) + radical(ROOJ) + radical(ROOJ)"""),
)

species(
    label = '[O]OC([O])[C]([O])OO(8906)',
    structure = SMILES('[O]OC([O])[C]([O])OO'),
    E0 = (134.811,'kJ/mol'),
    modes = [
        HarmonicOscillator(frequencies=([3615,1310,387.5,850,1000,360,370,350,1380,1390,370,380,2900,435,492.5,1135,1000,180,180,972.731],'cm^-1')),
        HinderedRotor(inertia=(0.1134,'amu*angstrom^2'), symmetry=1, barrier=(2.60728,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.113146,'amu*angstrom^2'), symmetry=1, barrier=(2.60145,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.113028,'amu*angstrom^2'), symmetry=1, barrier=(2.59873,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.113435,'amu*angstrom^2'), symmetry=1, barrier=(2.60809,'kJ/mol'), semiclassical=False),
    ],
    spinMultiplicity = 5,
    opticalIsomers = 1,
    molecularWeight = (122.034,'amu'),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[0.685054,0.0913309,-0.000191381,1.96187e-07,-7.37732e-11,16315.6,33.4259], Tmin=(100,'K'), Tmax=(855.674,'K')), NASAPolynomial(coeffs=[2.76868,0.0379896,-2.14407e-05,4.23438e-09,-2.92558e-13,17555.2,33.0255], Tmin=(855.674,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(134.811,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(216.176,'J/(mol*K)'), comment="""Thermo group additivity estimation: group(O2s-OsCs) + group(O2s-OsCs) + group(O2s-CsH) + group(O2s-CsH) + group(O2s-OsH) + group(O2s-OsH) + group(Cs-CsOsOsH) + group(Cs-CsOsOsH) + radical(CCOJ) + radical(CCOJ) + radical(ROOJ) + radical(Cs_P)"""),
)

species(
    label = '[O]O[C]([O])C(O)O[O](8907)',
    structure = SMILES('[O]O[C]([O])C(O)O[O]'),
    E0 = (61.1103,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
    molecularWeight = (122.034,'amu'),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[0.607375,0.0902088,-0.000182617,1.78466e-07,-6.43439e-11,7457.15,33.6051], Tmin=(100,'K'), Tmax=(870.231,'K')), NASAPolynomial(coeffs=[5.80398,0.0300313,-1.63353e-05,3.15924e-09,-2.14721e-13,7926.87,17.1523], Tmin=(870.231,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(61.1103,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(216.176,'J/(mol*K)'), comment="""Thermo group additivity estimation: group(O2s-OsCs) + group(O2s-OsCs) + group(O2s-CsH) + group(O2s-CsH) + group(O2s-OsH) + group(O2s-OsH) + group(Cs-CsOsOsH) + group(Cs-CsOsOsH) + radical(ROOJ) + radical(CCOJ) + radical(Cs_P) + radical(ROOJ)"""),
)

species(
    label = '[O]O[C]([O])C([O])OO(8908)',
    structure = SMILES('[O]O[C]([O])C([O])OO'),
    E0 = (134.811,'kJ/mol'),
    modes = [
        HarmonicOscillator(frequencies=([3615,1310,387.5,850,1000,360,370,350,1380,1390,370,380,2900,435,492.5,1135,1000,180,180,972.724],'cm^-1')),
        HinderedRotor(inertia=(0.113439,'amu*angstrom^2'), symmetry=1, barrier=(2.60818,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.113331,'amu*angstrom^2'), symmetry=1, barrier=(2.6057,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.113104,'amu*angstrom^2'), symmetry=1, barrier=(2.60049,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.113128,'amu*angstrom^2'), symmetry=1, barrier=(2.60102,'kJ/mol'), semiclassical=False),
    ],
    spinMultiplicity = 5,
    opticalIsomers = 1,
    molecularWeight = (122.034,'amu'),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[0.685193,0.0913291,-0.000191374,1.96176e-07,-7.37681e-11,16315.6,33.4254], Tmin=(100,'K'), Tmax=(855.69,'K')), NASAPolynomial(coeffs=[2.76837,0.0379901,-2.1441e-05,4.23446e-09,-2.92565e-13,17555.4,33.0272], Tmin=(855.69,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(134.811,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(216.176,'J/(mol*K)'), comment="""Thermo group additivity estimation: group(O2s-OsCs) + group(O2s-OsCs) + group(O2s-CsH) + group(O2s-CsH) + group(O2s-OsH) + group(O2s-OsH) + group(Cs-CsOsOsH) + group(Cs-CsOsOsH) + radical(Cs_P) + radical(CCOJ) + radical(CCOJ) + radical(ROOJ)"""),
)

species(
    label = '[O][CH]C([O])O[O](8416)',
    structure = SMILES('[O][CH]C([O])O[O]'),
    E0 = (238.588,'kJ/mol'),
    modes = [
        HarmonicOscillator(frequencies=([492.5,1135,1000,1380,1390,370,380,2900,435,3025,407.5,1350,352.5,180,1916.97,1924.25],'cm^-1')),
        HinderedRotor(inertia=(0.28295,'amu*angstrom^2'), symmetry=1, barrier=(6.50557,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.282123,'amu*angstrom^2'), symmetry=1, barrier=(6.48656,'kJ/mol'), semiclassical=False),
    ],
    spinMultiplicity = 5,
    opticalIsomers = 1,
    molecularWeight = (90.0349,'amu'),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[1.74873,0.0671988,-0.000152593,1.6348e-07,-6.22884e-11,28759.7,24.5506], Tmin=(100,'K'), Tmax=(875.926,'K')), NASAPolynomial(coeffs=[-0.266873,0.0323964,-1.76337e-05,3.40563e-09,-2.30945e-13,30801,43.6444], Tmin=(875.926,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(238.588,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(174.604,'J/(mol*K)'), comment="""Thermo group additivity estimation: group(O2s-OsCs) + group(O2s-CsH) + group(O2s-CsH) + group(O2s-OsH) + group(Cs-CsOsOsH) + group(Cs-CsOsHH) + radical(ROOJ) + radical(CCOJ) + radical(CCOJ) + radical(CCsJOH)"""),
)

species(
    label = '[O]O[C]([O])C([O])O[O](8909)',
    structure = SMILES('[O]O[C]([O])C([O])O[O]'),
    E0 = (286.815,'kJ/mol'),
    modes = [
        HarmonicOscillator(frequencies=([470,515,1100,1170,900,1100,1380,1390,370,380,2900,435,360,370,350,180,2012.64,2013.45],'cm^-1')),
        HinderedRotor(inertia=(0.245045,'amu*angstrom^2'), symmetry=1, barrier=(5.63408,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.245151,'amu*angstrom^2'), symmetry=1, barrier=(5.6365,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.245026,'amu*angstrom^2'), symmetry=1, barrier=(5.63363,'kJ/mol'), semiclassical=False),
    ],
    spinMultiplicity = 6,
    opticalIsomers = 1,
    molecularWeight = (121.026,'amu'),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[0.951555,0.0876908,-0.000195796,2.04848e-07,-7.70821e-11,34586,33.2292], Tmin=(100,'K'), Tmax=(870.409,'K')), NASAPolynomial(coeffs=[1.1344,0.0364104,-2.04979e-05,4.00531e-09,-2.73413e-13,36464.8,43.3484], Tmin=(870.409,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(286.815,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(195.39,'J/(mol*K)'), comment="""Thermo group additivity estimation: group(O2s-OsCs) + group(O2s-OsCs) + group(O2s-CsH) + group(O2s-CsH) + group(O2s-OsH) + group(O2s-OsH) + group(Cs-CsOsOsH) + group(Cs-CsOsOsH) + radical(CCOJ) + radical(Cs_P) + radical(ROOJ) + radical(CCOJ) + radical(ROOJ)"""),
)

species(
    label = '[O]OC([O])C(=O)OO(8910)',
    structure = SMILES('[O]OC([O])C(=O)OO'),
    E0 = (-338.416,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
    molecularWeight = (122.034,'amu'),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[1.10345,0.0622386,-6.82217e-05,3.67837e-08,-7.80517e-12,-40596.4,30.127], Tmin=(100,'K'), Tmax=(1146.63,'K')), NASAPolynomial(coeffs=[14.4133,0.0158067,-7.47951e-06,1.46683e-09,-1.04926e-13,-43648.7,-35.9064], Tmin=(1146.63,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(-338.416,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(216.176,'J/(mol*K)'), comment="""Thermo group additivity estimation: group(O2s-OsCs) + group(O2s-O2s(Cds-O2d)) + group(O2s-CsH) + group(O2s-OsH) + group(O2s-OsH) + group(Cs-CsOsOsH) + group(Cds-OdCsOs) + radical(C=OCOJ) + radical(ROOJ)"""),
)

species(
    label = '[O]OC(=O)C(O)O[O](8911)',
    structure = SMILES('[O]OC(=O)C(O)O[O]'),
    E0 = (-387.76,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
    molecularWeight = (122.034,'amu'),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[1.11162,0.0617588,-6.98586e-05,3.90545e-08,-8.55119e-12,-46531.2,30.0755], Tmin=(100,'K'), Tmax=(1116.79,'K')), NASAPolynomial(coeffs=[14.2772,0.0146046,-6.52528e-06,1.24844e-09,-8.82343e-14,-49471.9,-34.8953], Tmin=(1116.79,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(-387.76,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(216.176,'J/(mol*K)'), comment="""Thermo group additivity estimation: group(O2s-OsCs) + group(O2s-CsH) + group(O2s-O2s(Cds-O2d)) + group(O2s-OsH) + group(O2s-OsH) + group(Cs-CsOsOsH) + group(Cds-OdCsOs) + radical(C(=O)OOJ) + radical(ROOJ)"""),
)

species(
    label = '[O]OC(=O)C([O])OO(8912)',
    structure = SMILES('[O]OC(=O)C([O])OO'),
    E0 = (-296.032,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
    molecularWeight = (122.034,'amu'),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[1.10345,0.0622386,-6.82217e-05,3.67837e-08,-7.80517e-12,-35498.8,30.127], Tmin=(100,'K'), Tmax=(1146.63,'K')), NASAPolynomial(coeffs=[14.4133,0.0158067,-7.47951e-06,1.46683e-09,-1.04926e-13,-38551.1,-35.9064], Tmin=(1146.63,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(-296.032,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(216.176,'J/(mol*K)'), comment="""Thermo group additivity estimation: group(O2s-OsCs) + group(O2s-O2s(Cds-O2d)) + group(O2s-CsH) + group(O2s-OsH) + group(O2s-OsH) + group(Cs-CsOsOsH) + group(Cds-OdCsOs) + radical(C=OCOJ) + radical(C(=O)OOJ)"""),
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
    label = '[O]OC([O])C1OO1(8913)',
    structure = SMILES('[O]OC([O])C1OO1'),
    E0 = (26.5653,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
    molecularWeight = (106.034,'amu'),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[1.70695,0.0499821,-5.68543e-05,3.37979e-08,-7.97462e-12,3278.07,25.0036], Tmin=(100,'K'), Tmax=(1034.53,'K')), NASAPolynomial(coeffs=[10.7543,0.0150007,-6.13314e-06,1.11234e-09,-7.59184e-14,1406.14,-18.9515], Tmin=(1034.53,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(26.5653,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(199.547,'J/(mol*K)'), comment="""Thermo group additivity estimation: group(O2s-OsCs) + group(O2s-OsCs) + group(O2s-OsCs) + group(O2s-CsH) + group(O2s-OsH) + group(Cs-CsOsOsH) + group(Cs-CsOsOsH) + ring(dioxirane) + radical(CCOJ) + radical(ROOJ)"""),
)

species(
    label = '[O]OC1OOC1[O](8914)',
    structure = SMILES('[O]OC1OOC1[O]'),
    E0 = (39.0031,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
    molecularWeight = (106.034,'amu'),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[1.72638,0.0485218,-5.60743e-05,3.57194e-08,-8.9579e-12,4774.12,22.8914], Tmin=(100,'K'), Tmax=(1053.76,'K')), NASAPolynomial(coeffs=[9.70952,0.015361,-4.80323e-06,7.09265e-10,-4.13956e-14,3250.3,-15.2879], Tmin=(1053.76,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(39.0031,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(203.705,'J/(mol*K)'), comment="""Thermo group additivity estimation: group(O2s-OsCs) + group(O2s-OsCs) + group(O2s-OsCs) + group(O2s-CsH) + group(O2s-OsH) + group(Cs-CsOsOsH) + group(Cs-CsOsOsH) + ring(12dioxetane) + radical(ROOJ) + radical(CCOJ)"""),
)

species(
    label = '[O]C1OOOC1[O](8915)',
    structure = SMILES('[O]C1OOOC1[O]'),
    E0 = (1.60653,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
    molecularWeight = (106.034,'amu'),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[2.04566,0.0413364,-3.87265e-05,2.11454e-08,-4.76062e-12,264.984,20.8464], Tmin=(100,'K'), Tmax=(1068.94,'K')), NASAPolynomial(coeffs=[8.16933,0.0184211,-6.56986e-06,1.08994e-09,-7.00469e-14,-1044.16,-9.10493], Tmin=(1068.94,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(1.60653,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(207.862,'J/(mol*K)'), comment="""Thermo group additivity estimation: group(O2s-OsCs) + group(O2s-OsCs) + group(O2s-CsH) + group(O2s-CsH) + group(O2s-OsOs) + group(Cs-CsOsOsH) + group(Cs-CsOsOsH) + ring(123trioxolane) + radical(CCOJ) + radical(CCOJ)"""),
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
    label = '[O][CH]C(=O)O[O](8916)',
    structure = SMILES('[O][CH]C(=O)O[O]'),
    E0 = (-14.4008,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
    molecularWeight = (89.0269,'amu'),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[2.15356,0.0354757,-2.81459e-05,6.05093e-09,1.58788e-12,-1660.92,21.4889], Tmin=(100,'K'), Tmax=(972.295,'K')), NASAPolynomial(coeffs=[11.6239,0.00774599,-2.69235e-06,4.78395e-10,-3.39169e-14,-4033.37,-26.6637], Tmin=(972.295,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(-14.4008,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(149.66,'J/(mol*K)'), comment="""Thermo group additivity estimation: group(O2s-O2s(Cds-O2d)) + group(O2s-CsH) + group(O2s-OsH) + group(Cs-(Cds-O2d)OsHH) + group(Cds-OdCsOs) + radical(C=OCOJ) + radical(C(=O)OOJ) + radical(OCJC=O)"""),
)

species(
    label = '[O]OC([O])C1OOO1(8917)',
    structure = SMILES('[O]OC([O])C1OOO1'),
    E0 = (67.7777,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
    molecularWeight = (122.034,'amu'),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[1.46101,0.047539,-3.24162e-05,7.59904e-09,4.03599e-14,8250.01,28.9979], Tmin=(100,'K'), Tmax=(1271.73,'K')), NASAPolynomial(coeffs=[14.9489,0.016049,-8.17011e-06,1.64916e-09,-1.18975e-13,3935.27,-42.7917], Tmin=(1271.73,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(67.7777,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(224.491,'J/(mol*K)'), comment="""Thermo group additivity estimation: group(O2s-OsCs) + group(O2s-OsCs) + group(O2s-OsCs) + group(O2s-CsH) + group(O2s-OsOs) + group(O2s-OsH) + group(Cs-CsOsOsH) + group(Cs-CsOsOsH) + ring(Cyclobutane) + radical(CCOJ) + radical(ROOJ)"""),
)

species(
    label = '[O]OC1OOC1O[O](8196)',
    structure = SMILES('[O]OC1OOC1O[O]'),
    E0 = (36.8077,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
    molecularWeight = (122.034,'amu'),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[0.873508,0.0655919,-8.339e-05,5.36909e-08,-1.31572e-11,4542.3,25.9161], Tmin=(100,'K'), Tmax=(1114.35,'K')), NASAPolynomial(coeffs=[13.9609,0.0121661,-2.79516e-06,2.81881e-10,-1.01126e-14,2025.85,-36.8438], Tmin=(1114.35,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(36.8077,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(224.491,'J/(mol*K)'), comment="""Thermo group additivity estimation: group(O2s-OsCs) + group(O2s-OsCs) + group(O2s-OsCs) + group(O2s-OsCs) + group(O2s-OsH) + group(O2s-OsH) + group(Cs-CsOsOsH) + group(Cs-CsOsOsH) + ring(12dioxetane) + radical(ROOJ) + radical(ROOJ)"""),
)

species(
    label = '[O]OC1OOOC1[O](8918)',
    structure = SMILES('[O]OC1OOOC1[O]'),
    E0 = (-0.588926,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
    molecularWeight = (122.034,'amu'),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[1.35616,0.0563422,-5.81259e-05,2.77303e-08,-3.54313e-12,26.2023,24.68], Tmin=(100,'K'), Tmax=(841.415,'K')), NASAPolynomial(coeffs=[12.1808,0.0156635,-4.82594e-06,7.27006e-10,-4.42373e-14,-2177.02,-27.9413], Tmin=(841.415,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(-0.588926,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(228.648,'J/(mol*K)'), comment="""Thermo group additivity estimation: group(O2s-OsCs) + group(O2s-OsCs) + group(O2s-OsCs) + group(O2s-CsH) + group(O2s-OsOs) + group(O2s-OsH) + group(Cs-CsOsOsH) + group(Cs-CsOsOsH) + ring(123trioxolane) + radical(CCOJ) + radical(ROOJ)"""),
)

species(
    label = '[O]C1OOOOC1[O](8919)',
    structure = SMILES('[O]C1OOOOC1[O]'),
    E0 = (-2.23009,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
    molecularWeight = (122.034,'amu'),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[2.27323,0.0241908,3.57272e-05,-5.65202e-08,2.00473e-11,-194.474,24.109], Tmin=(100,'K'), Tmax=(1097.98,'K')), NASAPolynomial(coeffs=[12.0907,0.0241801,-1.31049e-05,2.78818e-09,-2.09732e-13,-4505.58,-33.987], Tmin=(1097.98,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(-2.23009,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(232.805,'J/(mol*K)'), comment="""Thermo group additivity estimation: group(O2s-OsCs) + group(O2s-OsCs) + group(O2s-CsH) + group(O2s-CsH) + group(O2s-OsOs) + group(O2s-OsOs) + group(Cs-CsOsOsH) + group(Cs-CsOsOsH) + ring(Cyclohexane) + radical(CCOJ) + radical(CCOJ)"""),
)

species(
    label = '[O]OC([O])C([O])[O](8920)',
    structure = SMILES('[O]OC([O])C([O])[O]'),
    E0 = (83.7646,'kJ/mol'),
    modes = [
        HarmonicOscillator(frequencies=([492.5,1135,1000,1380,1383.33,1386.67,1390,370,373.333,376.667,380,2800,3000,430,440,180,1939.97,1946.98,1948.02],'cm^-1')),
        HinderedRotor(inertia=(0.22477,'amu*angstrom^2'), symmetry=1, barrier=(5.16789,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.225173,'amu*angstrom^2'), symmetry=1, barrier=(5.17717,'kJ/mol'), semiclassical=False),
    ],
    spinMultiplicity = 5,
    opticalIsomers = 1,
    molecularWeight = (106.034,'amu'),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[1.69365,0.0708311,-0.000161743,1.80159e-07,-7.11889e-11,10138.3,28.5231], Tmin=(100,'K'), Tmax=(860.268,'K')), NASAPolynomial(coeffs=[-2.97661,0.0419964,-2.33242e-05,4.58606e-09,-3.16137e-13,12812.3,61.2232], Tmin=(860.268,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(83.7646,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(199.547,'J/(mol*K)'), comment="""Thermo group additivity estimation: group(O2s-OsCs) + group(O2s-CsH) + group(O2s-CsH) + group(O2s-CsH) + group(O2s-OsH) + group(Cs-CsOsOsH) + group(Cs-CsOsOsH) + radical(CCOJ) + radical(CCOJ) + radical(CCOJ) + radical(ROOJ)"""),
)

species(
    label = '[O]O[CH]C([O])O[O](8921)',
    structure = SMILES('[O]O[CH]C([O])O[O]'),
    E0 = (244.677,'kJ/mol'),
    modes = [
        HarmonicOscillator(frequencies=([470,515,1100,1170,900,1100,1380,1390,370,380,2900,435,3025,407.5,1350,352.5,180,1582.3],'cm^-1')),
        HinderedRotor(inertia=(0.259176,'amu*angstrom^2'), symmetry=1, barrier=(5.95896,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.259848,'amu*angstrom^2'), symmetry=1, barrier=(5.97441,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.26172,'amu*angstrom^2'), symmetry=1, barrier=(6.01746,'kJ/mol'), semiclassical=False),
    ],
    spinMultiplicity = 5,
    opticalIsomers = 1,
    molecularWeight = (106.034,'amu'),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[1.20171,0.0778277,-0.000163721,1.66353e-07,-6.13943e-11,29513.1,27.7744], Tmin=(100,'K'), Tmax=(876.632,'K')), NASAPolynomial(coeffs=[2.71725,0.0314274,-1.6763e-05,3.21325e-09,-2.17084e-13,30764.6,29.3158], Tmin=(876.632,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(244.677,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(195.39,'J/(mol*K)'), comment="""Thermo group additivity estimation: group(O2s-OsCs) + group(O2s-OsCs) + group(O2s-CsH) + group(O2s-OsH) + group(O2s-OsH) + group(Cs-CsOsOsH) + group(Cs-CsOsHH) + radical(CCsJOOH) + radical(ROOJ) + radical(CCOJ) + radical(ROOJ)"""),
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
    E0 = (81.5691,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS2',
    E0 = (114.283,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS3',
    E0 = (90.5476,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS4',
    E0 = (81.5691,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS5',
    E0 = (239.605,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS6',
    E0 = (269.983,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS7',
    E0 = (156.85,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS8',
    E0 = (218.7,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS9',
    E0 = (229.961,'kJ/mol'),
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
    E0 = (498.62,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS12',
    E0 = (144.969,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS13',
    E0 = (144.969,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS14',
    E0 = (106.542,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS15',
    E0 = (81.5691,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS16',
    E0 = (271.099,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS17',
    E0 = (283.132,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS18',
    E0 = (244.829,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS19',
    E0 = (220.938,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS20',
    E0 = (89.8534,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS21',
    E0 = (89.8534,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS22',
    E0 = (88.6819,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS23',
    E0 = (90.4392,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS24',
    E0 = (331.965,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS25',
    E0 = (651.559,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

reaction(
    label = 'reaction1',
    reactants = ['[O]OC([O])C([O])O[O](8193)'],
    products = ['[O]OC=O(5472)', '[O]OC=O(5472)'],
    transitionState = 'TS1',
    kinetics = Arrhenius(A=(5e+12,'s^-1'), n=0, Ea=(0,'kJ/mol'), T0=(1,'K'), Tmin=(300,'K'), Tmax=(1500,'K'), comment="""Exact match found for rate rule [RJJ]
Euclidian distance = 0
family: 1,4_Linear_birad_scission"""),
)

reaction(
    label = 'reaction2',
    reactants = ['H(8)', '[O]OC(=O)C([O])O[O](8904)'],
    products = ['[O]OC([O])C([O])O[O](8193)'],
    transitionState = 'TS2',
    kinetics = Arrhenius(A=(1.83701,'m^3/(mol*s)'), n=1.71338, Ea=(46.5052,'kJ/mol'), T0=(1,'K'), comment="""Estimated using average of templates [CO_O;HJ] + [CO-NdNd_O;YJ] for rate rule [CO-NdNd_O;HJ]
Euclidian distance = 1.0
family: R_Addition_MultipleBond"""),
)

reaction(
    label = 'reaction3',
    reactants = ['[O]OC=O(5472)', '[O][CH]O[O](8201)'],
    products = ['[O]OC([O])C([O])O[O](8193)'],
    transitionState = 'TS3',
    kinetics = Arrhenius(A=(6.87291,'m^3/(mol*s)'), n=1.39198, Ea=(57.1551,'kJ/mol'), T0=(1,'K'), comment="""Estimated using average of templates [CO_O;CJ] + [CO-NdH_O;YJ] for rate rule [CO-NdH_O;CJ]
Euclidian distance = 1.0
family: R_Addition_MultipleBond"""),
)

reaction(
    label = 'reaction4',
    reactants = ['O2(2)', '[O]OC([O])C=O(8158)'],
    products = ['[O]OC([O])C([O])O[O](8193)'],
    transitionState = 'TS4',
    kinetics = Arrhenius(A=(0.0849,'cm^3/(mol*s)'), n=3.486, Ea=(160.764,'kJ/mol'), T0=(1,'K'), Tmin=(300,'K'), Tmax=(1600,'K'), comment="""Estimated using template [CO-CsH_O;OJ] for rate rule [CO-CsH_O;O2b]
Euclidian distance = 1.0
Multiplied by reaction path degeneracy 2.0
family: R_Addition_MultipleBond
Ea raised from 158.3 to 160.8 kJ/mol to match endothermicity of reaction."""),
)

reaction(
    label = 'reaction5',
    reactants = ['[O]OC([O])C([O])O[O](8193)'],
    products = ['[O]O[C](O)C([O])O[O](8905)'],
    transitionState = 'TS5',
    kinetics = Arrhenius(A=(2.86762e+07,'s^-1'), n=1.70481, Ea=(158.036,'kJ/mol'), T0=(1,'K'), comment="""Estimated using template [R2H_S;Y_rad_out;Cs_H_out_NonDe] for rate rule [R2H_S;O_rad_out;Cs_H_out_NDMustO]
Euclidian distance = 1.41421356237
Multiplied by reaction path degeneracy 2.0
family: intra_H_migration"""),
)

reaction(
    label = 'reaction6',
    reactants = ['[O]OC([O])[C]([O])OO(8906)'],
    products = ['[O]OC([O])C([O])O[O](8193)'],
    transitionState = 'TS6',
    kinetics = Arrhenius(A=(40813.3,'s^-1'), n=2.17068, Ea=(135.172,'kJ/mol'), T0=(1,'K'), comment="""Estimated using average of templates [R3H_SS;Y_rad_out;O_H_out] + [R3H_SS_O;Y_rad_out;XH_out] for rate rule [R3H_SS_O;Y_rad_out;O_H_out]
Euclidian distance = 1.0
family: intra_H_migration"""),
)

reaction(
    label = 'reaction7',
    reactants = ['[O]OC([O])C([O])O[O](8193)'],
    products = ['[O]O[C]([O])C(O)O[O](8907)'],
    transitionState = 'TS7',
    kinetics = Arrhenius(A=(223829,'s^-1'), n=2.27675, Ea=(75.2806,'kJ/mol'), T0=(1,'K'), comment="""Estimated using an average for rate rule [R3H_SS_Cs;O_rad_out;XH_out]
Euclidian distance = 0
Multiplied by reaction path degeneracy 2.0
family: intra_H_migration"""),
)

reaction(
    label = 'reaction8',
    reactants = ['[O]O[C]([O])C([O])OO(8908)'],
    products = ['[O]OC([O])C([O])O[O](8193)'],
    transitionState = 'TS8',
    kinetics = Arrhenius(A=(0.00504539,'s^-1'), n=3.83, Ea=(83.8892,'kJ/mol'), T0=(1,'K'), comment="""Estimated using an average for rate rule [R4H_SSS;Y_rad_out;O_H_out]
Euclidian distance = 0
family: intra_H_migration"""),
)

reaction(
    label = 'reaction9',
    reactants = ['O2(2)', '[O][CH]C([O])O[O](8416)'],
    products = ['[O]OC([O])C([O])O[O](8193)'],
    transitionState = 'TS9',
    kinetics = Arrhenius(A=(2.18266e+06,'m^3/(mol*s)'), n=0.193158, Ea=(0,'kJ/mol'), T0=(1,'K'), comment="""Estimated using an average for rate rule [Y_rad;O2_birad]
Euclidian distance = 0
Multiplied by reaction path degeneracy 2.0
family: R_Recombination
Ea raised from -25.1 to 0 kJ/mol."""),
)

reaction(
    label = 'reaction10',
    reactants = ['[O][CH]O[O](8201)', '[O][CH]O[O](8201)'],
    products = ['[O]OC([O])C([O])O[O](8193)'],
    transitionState = 'TS10',
    kinetics = Arrhenius(A=(9.89449e+06,'m^3/(mol*s)'), n=-0.126319, Ea=(0,'kJ/mol'), T0=(1,'K'), comment="""Estimated using an average for rate rule [Y_rad;Y_rad]
Euclidian distance = 0
family: R_Recombination
Ea raised from -15.6 to -15.6 kJ/mol.
Ea raised from -15.6 to 0 kJ/mol."""),
)

reaction(
    label = 'reaction11',
    reactants = ['H(8)', '[O]O[C]([O])C([O])O[O](8909)'],
    products = ['[O]OC([O])C([O])O[O](8193)'],
    transitionState = 'TS11',
    kinetics = Arrhenius(A=(4.34078e+06,'m^3/(mol*s)'), n=0.278577, Ea=(0,'kJ/mol'), T0=(1,'K'), comment="""Estimated using an average for rate rule [Y_rad;H_rad]
Euclidian distance = 0
family: R_Recombination
Ea raised from -1.4 to 0 kJ/mol."""),
)

reaction(
    label = 'reaction12',
    reactants = ['[O]OC([O])C([O])O[O](8193)'],
    products = ['[O]OC([O])C(=O)OO(8910)'],
    transitionState = 'TS12',
    kinetics = Arrhenius(A=(1.4874e+09,'s^-1'), n=1.045, Ea=(63.4002,'kJ/mol'), T0=(1,'K'), comment="""Estimated using an average for rate rule [R3radExo;Y_rad;XH_Rrad]
Euclidian distance = 0
Multiplied by reaction path degeneracy 2.0
family: Intra_Disproportionation"""),
)

reaction(
    label = 'reaction13',
    reactants = ['[O]OC([O])C([O])O[O](8193)'],
    products = ['[O]OC(=O)C(O)O[O](8911)'],
    transitionState = 'TS13',
    kinetics = Arrhenius(A=(1.4874e+09,'s^-1'), n=1.045, Ea=(63.4002,'kJ/mol'), T0=(1,'K'), comment="""Estimated using an average for rate rule [R3radExo;Y_rad;XH_Rrad]
Euclidian distance = 0
Multiplied by reaction path degeneracy 2.0
family: Intra_Disproportionation"""),
)

reaction(
    label = 'reaction14',
    reactants = ['[O]OC([O])C([O])O[O](8193)'],
    products = ['[O]OC(=O)C([O])OO(8912)'],
    transitionState = 'TS14',
    kinetics = Arrhenius(A=(1.02844e+09,'s^-1'), n=0.311, Ea=(24.9733,'kJ/mol'), T0=(1,'K'), comment="""Estimated using template [R4;Y_rad;XH_Rrad] for rate rule [R4radExo;Y_rad;XH_Rrad]
Euclidian distance = 1.0
Multiplied by reaction path degeneracy 2.0
family: Intra_Disproportionation"""),
)

reaction(
    label = 'reaction15',
    reactants = ['[O]OC([O])C([O])O[O](8193)'],
    products = ['O2(S)(5486)', '[O]OC([O])C=O(8158)'],
    transitionState = 'TS15',
    kinetics = Arrhenius(A=(1e+13,'s^-1'), n=0, Ea=(0,'kJ/mol'), T0=(1,'K'), Tmin=(300,'K'), Tmax=(1500,'K'), comment="""Exact match found for rate rule [RJJ]
Euclidian distance = 0
Multiplied by reaction path degeneracy 2.0
family: 1,4_Linear_birad_scission"""),
)

reaction(
    label = 'reaction16',
    reactants = ['[O]OC([O])C([O])O[O](8193)'],
    products = ['O(T)(63)', '[O]OC([O])C1OO1(8913)'],
    transitionState = 'TS16',
    kinetics = Arrhenius(A=(4.60322e+09,'s^-1'), n=1.08, Ea=(189.53,'kJ/mol'), T0=(1,'K'), comment="""Estimated using an average for rate rule [R2OO_S;Y_rad_intra;OOJ]
Euclidian distance = 0
Multiplied by reaction path degeneracy 2.0
family: Cyclic_Ether_Formation"""),
)

reaction(
    label = 'reaction17',
    reactants = ['[O]OC([O])C([O])O[O](8193)'],
    products = ['O(T)(63)', '[O]OC1OOC1[O](8914)'],
    transitionState = 'TS17',
    kinetics = Arrhenius(A=(6.22709e+11,'s^-1'), n=0, Ea=(201.563,'kJ/mol'), T0=(1,'K'), comment="""Estimated using template [R3OO_SS;Y_rad_intra;OO] for rate rule [R3OO_SS;Y_rad_intra;OOJ]
Euclidian distance = 1.0
Multiplied by reaction path degeneracy 2.0
family: Cyclic_Ether_Formation"""),
)

reaction(
    label = 'reaction18',
    reactants = ['[O]OC([O])C([O])O[O](8193)'],
    products = ['O(T)(63)', '[O]C1OOOC1[O](8915)'],
    transitionState = 'TS18',
    kinetics = Arrhenius(A=(1.82826e+10,'s^-1'), n=0.54, Ea=(163.26,'kJ/mol'), T0=(1,'K'), comment="""Estimated using average of templates [RnOO;Y_rad_intra;OOJ] + [R4OO;Y_rad_intra;OO] for rate rule [R4OO;Y_rad_intra;OOJ]
Euclidian distance = 1.0
Multiplied by reaction path degeneracy 2.0
family: Cyclic_Ether_Formation"""),
)

reaction(
    label = 'reaction19',
    reactants = ['[O]OC([O])C([O])O[O](8193)'],
    products = ['[O]O(16)', '[O][CH]C(=O)O[O](8916)'],
    transitionState = 'TS19',
    kinetics = Arrhenius(A=(1.91635e+11,'s^-1'), n=0.573333, Ea=(139.369,'kJ/mol'), T0=(1,'K'), comment="""Estimated using an average for rate rule [R2OO]
Euclidian distance = 0
Multiplied by reaction path degeneracy 2.0
family: HO2_Elimination_from_PeroxyRadical"""),
)

reaction(
    label = 'reaction20',
    reactants = ['[O]OC([O])C([O])O[O](8193)'],
    products = ['[O]OC([O])C1OOO1(8917)'],
    transitionState = 'TS20',
    kinetics = Arrhenius(A=(3.24e+12,'s^-1'), n=-0.305, Ea=(8.28432,'kJ/mol'), T0=(1,'K'), comment="""Estimated using template [R4_SSS;Y_rad_out;Ypri_rad_out] for rate rule [R4_SSS;O_rad;Opri_rad]
Euclidian distance = 1.41421356237
Multiplied by reaction path degeneracy 2.0
family: Birad_recombination"""),
)

reaction(
    label = 'reaction21',
    reactants = ['[O]OC([O])C([O])O[O](8193)'],
    products = ['[O]OC1OOC1O[O](8196)'],
    transitionState = 'TS21',
    kinetics = Arrhenius(A=(1.62e+12,'s^-1'), n=-0.305, Ea=(8.28432,'kJ/mol'), T0=(1,'K'), comment="""Estimated using template [R4_SSS;Y_rad_out;Ypri_rad_out] for rate rule [R4_SSS;O_rad;Opri_rad]
Euclidian distance = 1.41421356237
family: Birad_recombination"""),
)

reaction(
    label = 'reaction22',
    reactants = ['[O]OC([O])C([O])O[O](8193)'],
    products = ['[O]OC1OOOC1[O](8918)'],
    transitionState = 'TS22',
    kinetics = Arrhenius(A=(1.552e+10,'s^-1'), n=0.311, Ea=(7.1128,'kJ/mol'), T0=(1,'K'), comment="""Estimated using template [R5_SSSS;Y_rad_out;Ypri_rad_out] for rate rule [R5_SSSS;O_rad;Opri_rad]
Euclidian distance = 1.41421356237
Multiplied by reaction path degeneracy 2.0
family: Birad_recombination"""),
)

reaction(
    label = 'reaction23',
    reactants = ['[O]OC([O])C([O])O[O](8193)'],
    products = ['[O]C1OOOOC1[O](8919)'],
    transitionState = 'TS23',
    kinetics = Arrhenius(A=(3.21e+10,'s^-1'), n=0.137, Ea=(8.87008,'kJ/mol'), T0=(1,'K'), comment="""Estimated using template [R6_SSSSS;Y_rad_out;Ypri_rad_out] for rate rule [R6_SSSSS;O_rad;Opri_rad]
Euclidian distance = 1.41421356237
family: Birad_recombination"""),
)

reaction(
    label = 'reaction24',
    reactants = ['O(T)(63)', '[O]OC([O])C([O])[O](8920)'],
    products = ['[O]OC([O])C([O])O[O](8193)'],
    transitionState = 'TS24',
    kinetics = Arrhenius(A=(87.1677,'m^3/(mol*s)'), n=1.88017, Ea=(5.1666,'kJ/mol'), T0=(1,'K'), Tmin=(303.03,'K'), Tmax=(2000,'K'), comment="""From training reaction 3 used for O_rad/NonDe;O_birad
Exact match found for rate rule [O_rad/NonDe;O_birad]
Euclidian distance = 0
Multiplied by reaction path degeneracy 2.0
family: Birad_R_Recombination"""),
)

reaction(
    label = 'reaction25',
    reactants = ['O(T)(63)', '[O]O[CH]C([O])O[O](8921)'],
    products = ['[O]OC([O])C([O])O[O](8193)'],
    transitionState = 'TS25',
    kinetics = Arrhenius(A=(93609.6,'m^3/(mol*s)'), n=1.13083, Ea=(163.847,'kJ/mol'), T0=(1,'K'), Tmin=(303.03,'K'), Tmax=(2000,'K'), comment="""Estimated using template [Y_rad;O_birad] for rate rule [C_rad/H/CsO;O_birad]
Euclidian distance = 4.0
family: Birad_R_Recombination"""),
)

network(
    label = '2307',
    isomers = [
        '[O]OC([O])C([O])O[O](8193)',
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
    label = '2307',
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

