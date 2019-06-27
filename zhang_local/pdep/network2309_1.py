species(
    label = '[O]O[CH]OO[CH]O[O](8195)',
    structure = SMILES('[O]O[CH]OO[CH]O[O]'),
    E0 = (293.432,'kJ/mol'),
    modes = [
        HarmonicOscillator(frequencies=([470,515,1100,1170,900,1100,350,500,795,815,3000,3050,390,425,1340,1360,335,370,180],'cm^-1')),
        HinderedRotor(inertia=(0.531278,'amu*angstrom^2'), symmetry=1, barrier=(12.2151,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.530932,'amu*angstrom^2'), symmetry=1, barrier=(12.2072,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(1.39908,'amu*angstrom^2'), symmetry=1, barrier=(32.1675,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(2.58241,'amu*angstrom^2'), symmetry=1, barrier=(59.3748,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.531548,'amu*angstrom^2'), symmetry=1, barrier=(12.2213,'kJ/mol'), semiclassical=False),
    ],
    spinMultiplicity = 5,
    opticalIsomers = 1,
    molecularWeight = (122.034,'amu'),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[-0.321025,0.106089,-0.000198394,1.74124e-07,-5.73022e-11,35436.9,30.9566], Tmin=(100,'K'), Tmax=(886.4,'K')), NASAPolynomial(coeffs=[13.7696,0.018576,-9.80999e-06,1.83618e-09,-1.21074e-13,33878.9,-30.0209], Tmin=(886.4,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(293.432,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(212.019,'J/(mol*K)'), comment="""Thermo group additivity estimation: group(O2s-OsCs) + group(O2s-OsCs) + group(O2s-OsCs) + group(O2s-OsCs) + group(O2s-OsH) + group(O2s-OsH) + group(Cs-OsOsHH) + group(Cs-OsOsHH) + radical(OCJO) + radical(OCJO) + radical(ROOJ) + radical(ROOJ)"""),
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
    label = '[O]O[CH]OOC=O(8882)',
    structure = SMILES('[O]O[CH]OOC=O'),
    E0 = (-169.447,'kJ/mol'),
    modes = [
        HarmonicOscillator(frequencies=([492.5,1135,1000,2782.5,750,1395,475,1775,1000,350,500,795,815,3025,407.5,1350,352.5],'cm^-1')),
        HinderedRotor(inertia=(1.66019,'amu*angstrom^2'), symmetry=1, barrier=(38.1709,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(1.66021,'amu*angstrom^2'), symmetry=1, barrier=(38.1714,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(1.66145,'amu*angstrom^2'), symmetry=1, barrier=(38.2,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(1.66297,'amu*angstrom^2'), symmetry=1, barrier=(38.2349,'kJ/mol'), semiclassical=False),
    ],
    spinMultiplicity = 3,
    opticalIsomers = 1,
    molecularWeight = (106.034,'amu'),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[0.994876,0.0579931,-5.6794e-05,2.33377e-08,-2.84291e-12,-20264.5,25.5526], Tmin=(100,'K'), Tmax=(1061.36,'K')), NASAPolynomial(coeffs=[17.3996,0.00829116,-3.68496e-06,7.40794e-10,-5.52988e-14,-24429.6,-57.7845], Tmin=(1061.36,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(-169.447,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(191.233,'J/(mol*K)'), comment="""Thermo group additivity estimation: group(O2s-OsCs) + group(O2s-OsCs) + group(O2s-O2s(Cds-O2d)) + group(O2s-OsH) + group(Cs-OsOsHH) + group(Cds-OdOsH) + radical(OCJO) + radical(ROOJ)"""),
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
    label = '[O]O[CH]OOC1OO1(8883)',
    structure = SMILES('[O]O[CH]OOC1OO1'),
    E0 = (85.4496,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
    molecularWeight = (122.034,'amu'),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[0.700708,0.0766846,-0.000114856,8.78045e-08,-2.63816e-11,10392.2,27.3281], Tmin=(100,'K'), Tmax=(819.468,'K')), NASAPolynomial(coeffs=[12.544,0.0188748,-9.03733e-06,1.71685e-09,-1.1829e-13,8451.18,-27.451], Tmin=(819.468,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(85.4496,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(216.176,'J/(mol*K)'), comment="""Thermo group additivity estimation: group(O2s-OsCs) + group(O2s-OsCs) + group(O2s-OsCs) + group(O2s-OsCs) + group(O2s-OsCs) + group(O2s-OsH) + group(Cs-OsOsOsH) + group(Cs-OsOsHH) + ring(dioxirane) + radical(ROOJ) + radical(OCJO)"""),
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
    label = '[O]O[CH]O[O](8054)',
    structure = SMILES('[O]O[CH]O[O]'),
    E0 = (225,'kJ/mol'),
    modes = [
        HarmonicOscillator(frequencies=([470,515,1100,1170,900,1100,3025,407.5,1350,352.5],'cm^-1')),
        HinderedRotor(inertia=(0.512555,'amu*angstrom^2'), symmetry=1, barrier=(11.7847,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.513195,'amu*angstrom^2'), symmetry=1, barrier=(11.7994,'kJ/mol'), semiclassical=False),
    ],
    spinMultiplicity = 4,
    opticalIsomers = 1,
    molecularWeight = (77.0162,'amu'),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[1.97108,0.0511084,-9.91678e-05,8.86788e-08,-2.91459e-11,27128.2,17.3189], Tmin=(100,'K'), Tmax=(920.661,'K')), NASAPolynomial(coeffs=[7.87783,0.00915139,-4.26167e-06,7.32392e-10,-4.46318e-14,26731.1,-6.93917], Tmin=(920.661,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(225,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(124.717,'J/(mol*K)'), comment="""Thermo group additivity estimation: group(O2s-OsCs) + group(O2s-OsCs) + group(O2s-OsH) + group(O2s-OsH) + group(Cs-OsOsHH) + radical(ROOJ) + radical(ROOJ) + radical(OCJO)"""),
)

species(
    label = '[O][CH]OO[CH]O[O](8884)',
    structure = SMILES('[O][CH]OO[CH]O[O]'),
    E0 = (297.321,'kJ/mol'),
    modes = [
        HarmonicOscillator(frequencies=([492.5,1135,1000,350,500,795,815,3000,3050,390,425,1340,1360,335,370,180,1378.69],'cm^-1')),
        HinderedRotor(inertia=(0.296884,'amu*angstrom^2'), symmetry=1, barrier=(6.82595,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.296061,'amu*angstrom^2'), symmetry=1, barrier=(6.80703,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(2.28465,'amu*angstrom^2'), symmetry=1, barrier=(52.5286,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(2.28489,'amu*angstrom^2'), symmetry=1, barrier=(52.534,'kJ/mol'), semiclassical=False),
    ],
    spinMultiplicity = 5,
    opticalIsomers = 1,
    molecularWeight = (106.034,'amu'),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[0.676084,0.0895528,-0.000185787,1.84284e-07,-6.71865e-11,35863.4,26.3298], Tmin=(100,'K'), Tmax=(866.889,'K')), NASAPolynomial(coeffs=[5.10733,0.0305971,-1.71417e-05,3.34636e-09,-2.28608e-13,36542,13.9303], Tmin=(866.889,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(297.321,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(191.233,'J/(mol*K)'), comment="""Thermo group additivity estimation: group(O2s-OsCs) + group(O2s-OsCs) + group(O2s-OsCs) + group(O2s-CsH) + group(O2s-OsH) + group(Cs-OsOsHH) + group(Cs-OsOsHH) + radical(OCJO) + radical(OCOJ) + radical(OCJO) + radical(ROOJ)"""),
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
    label = '[CH]OO[CH]O[O](8389)',
    structure = SMILES('[CH]OO[CH]O[O]'),
    E0 = (533.908,'kJ/mol'),
    modes = [
        HarmonicOscillator(frequencies=([492.5,1135,1000,350,500,795,815,3025,407.5,1350,352.5,345.394,345.395,1639.04],'cm^-1')),
        HinderedRotor(inertia=(0.102518,'amu*angstrom^2'), symmetry=1, barrier=(8.67879,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.102518,'amu*angstrom^2'), symmetry=1, barrier=(8.67879,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.102518,'amu*angstrom^2'), symmetry=1, barrier=(8.67879,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.465624,'amu*angstrom^2'), symmetry=1, barrier=(39.4177,'kJ/mol'), semiclassical=False),
    ],
    spinMultiplicity = 5,
    opticalIsomers = 1,
    molecularWeight = (90.0349,'amu'),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[0.983951,0.0734274,-0.000133554,1.16993e-07,-3.90182e-11,64316.1,22.8819], Tmin=(100,'K'), Tmax=(854.454,'K')), NASAPolynomial(coeffs=[10.7545,0.0148503,-8.18438e-06,1.59276e-09,-1.08934e-13,63115,-19.9759], Tmin=(854.454,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(533.908,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(166.289,'J/(mol*K)'), comment="""Thermo group additivity estimation: group(O2s-OsCs) + group(O2s-OsCs) + group(O2s-OsCs) + group(O2s-OsH) + group(Cs-OsOsHH) + group(Cs-OsHHH) + radical(ROOJ) + radical(CH2_triplet) + radical(OCJO)"""),
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
    label = '[O]O[C]OO[CH]O[O](8885)',
    structure = SMILES('[O]O[C]OO[CH]O[O]'),
    E0 = (565.461,'kJ/mol'),
    modes = [
        HarmonicOscillator(frequencies=([470,515,1100,1170,900,1100,350,500,795,815,3025,407.5,1350,352.5,180,180],'cm^-1')),
        HinderedRotor(inertia=(1.02258,'amu*angstrom^2'), symmetry=1, barrier=(23.5112,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(1.02143,'amu*angstrom^2'), symmetry=1, barrier=(23.4846,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(1.02389,'amu*angstrom^2'), symmetry=1, barrier=(23.5412,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(1.02414,'amu*angstrom^2'), symmetry=1, barrier=(23.5471,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(1.02461,'amu*angstrom^2'), symmetry=1, barrier=(23.5577,'kJ/mol'), semiclassical=False),
    ],
    spinMultiplicity = 6,
    opticalIsomers = 1,
    molecularWeight = (121.026,'amu'),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[-0.550141,0.0964599,-0.000145691,9.69256e-08,-2.40637e-11,68176.3,30.72], Tmin=(100,'K'), Tmax=(1009.41,'K')), NASAPolynomial(coeffs=[24.4552,-0.00263076,1.56214e-06,-3.29275e-10,2.37362e-14,63128.2,-90.1503], Tmin=(1009.41,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(565.461,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(187.076,'J/(mol*K)'), comment="""Thermo group additivity estimation: group(O2s-OsCs) + group(O2s-OsCs) + group(O2s-OsCs) + group(O2s-OsCs) + group(O2s-OsH) + group(O2s-OsH) + group(Cs-OsOsHH) + group(Cs-OsOsHH) + radical(OCJO) + radical(ROOJ) + radical(CH2_triplet) + radical(ROOJ)"""),
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
    E0 = (293.432,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS2',
    E0 = (293.856,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS3',
    E0 = (457.776,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS4',
    E0 = (301.293,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS5',
    E0 = (301.717,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS6',
    E0 = (695.641,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS7',
    E0 = (704.202,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS8',
    E0 = (559.596,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS9',
    E0 = (777.265,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

reaction(
    label = 'reaction1',
    reactants = ['[O]O[CH]OO[CH]O[O](8195)'],
    products = ['[O]OC=O(5472)', '[O]OC=O(5472)'],
    transitionState = 'TS1',
    kinetics = Arrhenius(A=(5e+12,'s^-1'), n=0, Ea=(0,'kJ/mol'), T0=(1,'K'), Tmin=(300,'K'), Tmax=(1500,'K'), comment="""Exact match found for rate rule [RJJ]
Euclidian distance = 0
family: 1,4_Linear_birad_scission"""),
)

reaction(
    label = 'reaction2',
    reactants = ['O(T)(63)', '[O]O[CH]OOC=O(8882)'],
    products = ['[O]O[CH]OO[CH]O[O](8195)'],
    transitionState = 'TS2',
    kinetics = Arrhenius(A=(4000,'m^3/(mol*s)'), n=1.39, Ea=(220.269,'kJ/mol'), T0=(1,'K'), comment="""Estimated using template [Od_CO-NdH;YJ] for rate rule [Od_CO-NdH;O_atom_triplet]
Euclidian distance = 2.0
family: R_Addition_MultipleBond"""),
)

reaction(
    label = 'reaction3',
    reactants = ['[O][CH]O[O](8201)', '[O][CH]O[O](8201)'],
    products = ['[O]O[CH]OO[CH]O[O](8195)'],
    transitionState = 'TS3',
    kinetics = Arrhenius(A=(9.89449e+06,'m^3/(mol*s)'), n=-0.126319, Ea=(0,'kJ/mol'), T0=(1,'K'), comment="""Estimated using an average for rate rule [Y_rad;Y_rad]
Euclidian distance = 0
family: R_Recombination
Ea raised from -15.6 to -15.6 kJ/mol.
Ea raised from -15.6 to 0 kJ/mol."""),
)

reaction(
    label = 'reaction4',
    reactants = ['[O]O[CH]OO[CH]O[O](8195)'],
    products = ['[O]O[CH]OOC1OO1(8883)'],
    transitionState = 'TS4',
    kinetics = Arrhenius(A=(3.23934e+11,'s^-1'), n=0.0247333, Ea=(7.86034,'kJ/mol'), T0=(1,'K'), comment="""Estimated using template [Rn;Y_rad_out;Cpri_rad_out_single] for rate rule [R3_SS;O_rad;Cpri_rad_out_H/NonDeO]
Euclidian distance = 3.0
Multiplied by reaction path degeneracy 2.0
family: Birad_recombination"""),
)

reaction(
    label = 'reaction5',
    reactants = ['[O]O[CH]OO[CH]O[O](8195)'],
    products = ['[O]OC1OOC1O[O](8196)'],
    transitionState = 'TS5',
    kinetics = Arrhenius(A=(1.62e+12,'s^-1'), n=-0.305, Ea=(8.28432,'kJ/mol'), T0=(1,'K'), comment="""Estimated using template [R4_SSS;C_rad_out_single;Cpri_rad_out_single] for rate rule [R4_SSS;C_rad_out_H/NonDeO;Cpri_rad_out_H/NonDeO]
Euclidian distance = 2.82842712475
family: Birad_recombination"""),
)

reaction(
    label = 'reaction6',
    reactants = ['[CH]O[O](2819)', '[O]O[CH]O[O](8054)'],
    products = ['[O]O[CH]OO[CH]O[O](8195)'],
    transitionState = 'TS6',
    kinetics = Arrhenius(A=(87.1677,'m^3/(mol*s)'), n=1.88017, Ea=(5.1666,'kJ/mol'), T0=(1,'K'), Tmin=(303.03,'K'), Tmax=(2000,'K'), comment="""Estimated using an average for rate rule [O_rad/NonDe;Birad]
Euclidian distance = 0
Multiplied by reaction path degeneracy 2.0
family: Birad_R_Recombination"""),
)

reaction(
    label = 'reaction7',
    reactants = ['O(T)(63)', '[O][CH]OO[CH]O[O](8884)'],
    products = ['[O]O[CH]OO[CH]O[O](8195)'],
    transitionState = 'TS7',
    kinetics = Arrhenius(A=(93609.6,'m^3/(mol*s)'), n=1.13083, Ea=(163.847,'kJ/mol'), T0=(1,'K'), Tmin=(303.03,'K'), Tmax=(2000,'K'), comment="""From training reaction 2 used for Y_rad;O_birad
Exact match found for rate rule [Y_rad;O_birad]
Euclidian distance = 0
family: Birad_R_Recombination"""),
)

reaction(
    label = 'reaction8',
    reactants = ['O2(2)', '[CH]OO[CH]O[O](8389)'],
    products = ['[O]O[CH]OO[CH]O[O](8195)'],
    transitionState = 'TS8',
    kinetics = Arrhenius(A=(2.29708e+06,'m^3/(mol*s)'), n=0.575199, Ea=(34.3157,'kJ/mol'), T0=(1,'K'), comment="""Estimated using an average for rate rule [Y_rad;Birad]
Euclidian distance = 0
Multiplied by reaction path degeneracy 2.0
family: Birad_R_Recombination"""),
)

reaction(
    label = 'reaction9',
    reactants = ['H(8)', '[O]O[C]OO[CH]O[O](8885)'],
    products = ['[O]O[CH]OO[CH]O[O](8195)'],
    transitionState = 'TS9',
    kinetics = Arrhenius(A=(1e+07,'m^3/(mol*s)'), n=0, Ea=(0,'kJ/mol'), T0=(1,'K'), comment="""Estimated using an average for rate rule [H_rad;Birad]
Euclidian distance = 0
family: Birad_R_Recombination"""),
)

network(
    label = '2309',
    isomers = [
        '[O]O[CH]OO[CH]O[O](8195)',
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
    label = '2309',
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

