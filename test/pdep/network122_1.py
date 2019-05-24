species(
    label = '[CH2][CH]CC([CH2])[CH2](214)',
    structure = SMILES('[CH2][CH]CC([CH2])[CH2]'),
    E0 = (612.298,'kJ/mol'),
    modes = [
        HarmonicOscillator(frequencies=([2750,2850,1437.5,1250,1305,750,350,1380,1390,370,380,2900,435,3000,3020,3040,3060,3080,3100,415,440,465,780,815,850,1435,1455,1475,900,1000,1100,3025,407.5,1350,352.5,1170.55,1171.88],'cm^-1')),
        HinderedRotor(inertia=(0.168703,'amu*angstrom^2'), symmetry=1, barrier=(3.87882,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.00398897,'amu*angstrom^2'), symmetry=1, barrier=(3.86095,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.168672,'amu*angstrom^2'), symmetry=1, barrier=(3.87809,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.171793,'amu*angstrom^2'), symmetry=1, barrier=(3.94985,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.00403782,'amu*angstrom^2'), symmetry=1, barrier=(3.94605,'kJ/mol'), semiclassical=False),
    ],
    spinMultiplicity = 5,
    opticalIsomers = 1,
    molecularWeight = (82.1436,'amu'),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[1.64318,0.0551015,-3.22854e-05,-1.62113e-09,1.00768e-11,73724.5,29.706], Tmin=(100,'K'), Tmax=(647.81,'K')), NASAPolynomial(coeffs=[5.92381,0.0380635,-1.45841e-05,2.54567e-09,-1.69303e-13,72972.8,9.39165], Tmin=(647.81,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(612.298,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(361.68,'J/(mol*K)'), comment="""Thermo group additivity estimation: group(Cs-CsCsCsH) + group(Cs-CsCsHH) + longDistanceInteraction_noncyclic(CsCs-ST) + group(Cs-CsCsHH) + group(Cs-CsHHH) + group(Cs-CsHHH) + group(Cs-CsHHH) + radical(Isobutyl) + radical(Isobutyl) + radical(RCCJC) + radical(RCCJ)"""),
)

species(
    label = '[CH2]C=C(66)',
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
    label = 'H(7)',
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
    label = '[CH2][CH]CC([CH2])=C(292)',
    structure = SMILES('[CH2][CH]CC([CH2])=C'),
    E0 = (471.835,'kJ/mol'),
    modes = [
        HarmonicOscillator(frequencies=([2750,2850,1437.5,1250,1305,750,350,2950,3100,1380,975,1025,1650,350,440,435,1725,3000,3033.33,3066.67,3100,415,465,780,850,1435,1475,900,1100,3025,407.5,1350,352.5,180,1119.79],'cm^-1')),
        HinderedRotor(inertia=(0.004173,'amu*angstrom^2'), symmetry=1, barrier=(3.71386,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.161098,'amu*angstrom^2'), symmetry=1, barrier=(3.70396,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.161509,'amu*angstrom^2'), symmetry=1, barrier=(3.71341,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.0642308,'amu*angstrom^2'), symmetry=1, barrier=(57.1255,'kJ/mol'), semiclassical=False),
    ],
    spinMultiplicity = 4,
    opticalIsomers = 1,
    molecularWeight = (81.1357,'amu'),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[1.42509,0.0531145,-3.31274e-05,1.07407e-08,-1.45658e-12,56843.8,26.433], Tmin=(100,'K'), Tmax=(1641.37,'K')), NASAPolynomial(coeffs=[10.9944,0.0297939,-1.18149e-05,2.08413e-09,-1.38062e-13,53702.6,-24.4748], Tmin=(1641.37,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(471.835,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(340.893,'J/(mol*K)'), comment="""Thermo group additivity estimation: group(Cs-CsCsHH) + group(Cs-(Cds-Cds)CsHH) + group(Cs-CsHHH) + group(Cs-(Cds-Cds)HHH) + group(Cds-CdsCsCs) + group(Cds-CdsHH) + radical(RCCJC) + radical(Allyl_P) + radical(RCCJ)"""),
)

species(
    label = '[CH2]C=CC([CH2])[CH2](293)',
    structure = SMILES('[CH2]C=CC([CH2])[CH2]'),
    E0 = (477.364,'kJ/mol'),
    modes = [
        HarmonicOscillator(frequencies=([1380,1390,370,380,2900,435,3000,3020,3040,3060,3080,3100,415,440,465,780,815,850,1435,1455,1475,900,1000,1100,2995,3025,975,1000,1300,1375,400,500,1630,1680,246.839],'cm^-1')),
        HinderedRotor(inertia=(0.00276675,'amu*angstrom^2'), symmetry=1, barrier=(0.119627,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(1.71069,'amu*angstrom^2'), symmetry=1, barrier=(73.9598,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.00276659,'amu*angstrom^2'), symmetry=1, barrier=(0.119627,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.0361335,'amu*angstrom^2'), symmetry=1, barrier=(73.9606,'kJ/mol'), semiclassical=False),
    ],
    spinMultiplicity = 4,
    opticalIsomers = 1,
    molecularWeight = (81.1357,'amu'),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[1.36913,0.0491307,-9.81018e-06,-2.3684e-08,1.37807e-11,57516.4,26.4779], Tmin=(100,'K'), Tmax=(924.202,'K')), NASAPolynomial(coeffs=[11.9922,0.0255949,-8.03368e-06,1.30778e-09,-8.67286e-14,54594.4,-29.1199], Tmin=(924.202,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(477.364,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(340.893,'J/(mol*K)'), comment="""Thermo group additivity estimation: group(Cs-(Cds-Cds)CsCsH) + group(Cs-CsHHH) + group(Cs-CsHHH) + group(Cs-(Cds-Cds)HHH) + group(Cds-CdsCsH) + group(Cds-CdsCsH) + radical(Isobutyl) + radical(Isobutyl) + radical(Allyl_P)"""),
)

species(
    label = '[CH2][CH]CC=C(188)',
    structure = SMILES('[CH2][CH]CC=C'),
    E0 = (359.391,'kJ/mol'),
    modes = [
        HarmonicOscillator(frequencies=([2750,2850,1437.5,1250,1305,750,350,2950,3100,1380,975,1025,1650,3010,987.5,1337.5,450,1655,3000,3100,440,815,1455,1000,3025,407.5,1350,352.5,2209.08,2209.39],'cm^-1')),
        HinderedRotor(inertia=(0.0586692,'amu*angstrom^2'), symmetry=1, barrier=(6.0365,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.0586378,'amu*angstrom^2'), symmetry=1, barrier=(6.03061,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.401657,'amu*angstrom^2'), symmetry=1, barrier=(41.4786,'kJ/mol'), semiclassical=False),
    ],
    spinMultiplicity = 3,
    opticalIsomers = 1,
    molecularWeight = (68.117,'amu'),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[2.13435,0.0369817,-1.51233e-05,1.33711e-09,3.7705e-13,43295.2,23.5599], Tmin=(100,'K'), Tmax=(1532.35,'K')), NASAPolynomial(coeffs=[8.68443,0.0259362,-1.02359e-05,1.7885e-09,-1.17141e-13,40577.2,-13.1548], Tmin=(1532.35,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(359.391,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(295.164,'J/(mol*K)'), comment="""Thermo group additivity estimation: group(Cs-CsCsHH) + group(Cs-(Cds-Cds)CsHH) + group(Cs-CsHHH) + group(Cds-CdsCsH) + group(Cds-CdsHH) + radical(RCCJC) + radical(RCCJ)"""),
)

species(
    label = 'CH2(T)(26)',
    structure = SMILES('[CH2]'),
    E0 = (381.37,'kJ/mol'),
    modes = [
        HarmonicOscillator(frequencies=([1066.91,2790.99,3622.37],'cm^-1')),
    ],
    spinMultiplicity = 3,
    opticalIsomers = 1,
    molecularWeight = (14.0266,'amu'),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[4.01192,-0.000154979,3.26298e-06,-2.40422e-09,5.69497e-13,45867.7,0.5332], Tmin=(100,'K'), Tmax=(1104.58,'K')), NASAPolynomial(coeffs=[3.14983,0.00296674,-9.76056e-07,1.54115e-10,-9.50338e-15,46058.1,4.77808], Tmin=(1104.58,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(381.37,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(58.2013,'J/(mol*K)'), label="""CH2(T)""", comment="""Thermo library: primaryThermoLibrary"""),
)

species(
    label = '[CH2][CH][CH2](219)',
    structure = SMILES('[CH2][CH][CH2]'),
    E0 = (484.293,'kJ/mol'),
    modes = [
        HarmonicOscillator(frequencies=([3025,407.5,1350,352.5,3000,3033.33,3066.67,3100,415,465,780,850,1435,1475,900,1100],'cm^-1')),
        HinderedRotor(inertia=(0.00132734,'amu*angstrom^2'), symmetry=1, barrier=(2.4107,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.00132734,'amu*angstrom^2'), symmetry=1, barrier=(2.4107,'kJ/mol'), semiclassical=False),
    ],
    spinMultiplicity = 4,
    opticalIsomers = 1,
    molecularWeight = (41.0718,'amu'),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[3.34321,0.013802,2.16438e-06,-5.76341e-09,1.61336e-12,58271,14.9549], Tmin=(100,'K'), Tmax=(1447.1,'K')), NASAPolynomial(coeffs=[4.39495,0.0167646,-6.99098e-06,1.25743e-09,-8.3812e-14,57352,7.36867], Tmin=(1447.1,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(484.293,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(174.604,'J/(mol*K)'), comment="""Thermo group additivity estimation: group(Cs-CsCsHH) + group(Cs-CsHHH) + group(Cs-CsHHH) + radical(CCJC) + radical(RCCJ) + radical(RCCJ)"""),
)

species(
    label = '[CH2][CH]C[C]([CH2])C(257)',
    structure = SMILES('[CH2][CH]C[C]([CH2])C'),
    E0 = (592.639,'kJ/mol'),
    modes = [
        HarmonicOscillator(frequencies=([3000,3033.33,3066.67,3100,415,465,780,850,1435,1475,900,1100,2750,2850,1437.5,1250,1305,750,350,2750,2800,2850,1350,1500,750,1050,1375,1000,3025,407.5,1350,352.5,360,370,350,1668.23,1668.63],'cm^-1')),
        HinderedRotor(inertia=(0.00254354,'amu*angstrom^2'), symmetry=1, barrier=(5.02347,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.218446,'amu*angstrom^2'), symmetry=1, barrier=(5.02251,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.00254252,'amu*angstrom^2'), symmetry=1, barrier=(5.02098,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.218389,'amu*angstrom^2'), symmetry=1, barrier=(5.02119,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.218232,'amu*angstrom^2'), symmetry=1, barrier=(5.01758,'kJ/mol'), semiclassical=False),
    ],
    spinMultiplicity = 5,
    opticalIsomers = 1,
    molecularWeight = (82.1436,'amu'),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[1.43297,0.0648349,-8.48135e-05,7.94924e-08,-2.97461e-11,71362.4,30.4654], Tmin=(100,'K'), Tmax=(864.667,'K')), NASAPolynomial(coeffs=[0.922137,0.0475416,-2.0714e-05,3.77992e-09,-2.53872e-13,72185.5,37.1047], Tmin=(864.667,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(592.639,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(361.68,'J/(mol*K)'), comment="""Thermo group additivity estimation: group(Cs-CsCsCsH) + group(Cs-CsCsHH) + longDistanceInteraction_noncyclic(CsCs-ST) + group(Cs-CsCsHH) + group(Cs-CsHHH) + group(Cs-CsHHH) + group(Cs-CsHHH) + radical(RCCJ) + radical(Tertalkyl) + radical(Isobutyl) + radical(RCCJC)"""),
)

species(
    label = '[CH2]C[CH]C([CH2])[CH2](294)',
    structure = SMILES('[CH2]C[CH]C([CH2])[CH2]'),
    E0 = (612.394,'kJ/mol'),
    modes = [
        HarmonicOscillator(frequencies=([2750,2850,1437.5,1250,1305,750,350,1380,1390,370,380,2900,435,3000,3020,3040,3060,3080,3100,415,440,465,780,815,850,1435,1455,1475,900,1000,1100,3025,407.5,1350,352.5,1474.84,1474.86],'cm^-1')),
        HinderedRotor(inertia=(0.00293216,'amu*angstrom^2'), symmetry=1, barrier=(4.52599,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.00293217,'amu*angstrom^2'), symmetry=1, barrier=(4.52599,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(2.83834,'amu*angstrom^2'), symmetry=1, barrier=(65.2589,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.196827,'amu*angstrom^2'), symmetry=1, barrier=(4.52543,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.196846,'amu*angstrom^2'), symmetry=1, barrier=(4.52587,'kJ/mol'), semiclassical=False),
    ],
    spinMultiplicity = 5,
    opticalIsomers = 1,
    molecularWeight = (82.1436,'amu'),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[1.34718,0.0577643,-4.24404e-05,1.89838e-08,-3.75199e-12,73750,30.5481], Tmin=(100,'K'), Tmax=(1153.12,'K')), NASAPolynomial(coeffs=[7.69479,0.0357454,-1.37978e-05,2.42436e-09,-1.61845e-13,72286.1,-0.979865], Tmin=(1153.12,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(612.394,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(361.68,'J/(mol*K)'), comment="""Thermo group additivity estimation: group(Cs-CsCsCsH) + group(Cs-CsCsHH) + longDistanceInteraction_noncyclic(CsCs-ST) + group(Cs-CsCsHH) + group(Cs-CsHHH) + group(Cs-CsHHH) + group(Cs-CsHHH) + radical(Cs_S) + radical(Isobutyl) + radical(Isobutyl) + radical(RCCJ)"""),
)

species(
    label = '[CH2][CH][CH]C([CH2])C(258)',
    structure = SMILES('[CH2][CH][CH]C([CH2])C'),
    E0 = (601.758,'kJ/mol'),
    modes = [
        HarmonicOscillator(frequencies=([3000,3050,390,425,1340,1360,335,370,1380,1390,370,380,2900,435,3000,3033.33,3066.67,3100,415,465,780,850,1435,1475,900,1100,2750,2800,2850,1350,1500,750,1050,1375,1000,1360.78,1362.84],'cm^-1')),
        HinderedRotor(inertia=(0.165884,'amu*angstrom^2'), symmetry=1, barrier=(3.81401,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.00288411,'amu*angstrom^2'), symmetry=1, barrier=(3.7978,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.165332,'amu*angstrom^2'), symmetry=1, barrier=(3.8013,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.00289531,'amu*angstrom^2'), symmetry=1, barrier=(3.8125,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.166177,'amu*angstrom^2'), symmetry=1, barrier=(3.82074,'kJ/mol'), semiclassical=False),
    ],
    spinMultiplicity = 5,
    opticalIsomers = 1,
    molecularWeight = (82.1436,'amu'),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[2.4336,0.041928,2.62665e-05,-1.26438e-07,1.09331e-10,72423.3,28.4918], Tmin=(100,'K'), Tmax=(437.156,'K')), NASAPolynomial(coeffs=[3.59069,0.0431409,-1.83848e-05,3.40343e-09,-2.34319e-13,72209.4,22.5773], Tmin=(437.156,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(601.758,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(361.68,'J/(mol*K)'), comment="""Thermo group additivity estimation: group(Cs-CsCsCsH) + group(Cs-CsCsHH) + longDistanceInteraction_noncyclic(CsCs-ST) + group(Cs-CsCsHH) + group(Cs-CsHHH) + group(Cs-CsHHH) + group(Cs-CsHHH) + radical(Cs_S) + radical(RCCJ) + radical(RCCJC) + radical(Isobutyl)"""),
)

species(
    label = '[CH2]CC[C]([CH2])[CH2](295)',
    structure = SMILES('[CH2]CC[C]([CH2])[CH2]'),
    E0 = (603.275,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
    molecularWeight = (82.1436,'amu'),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[1.19535,0.0669514,-7.98698e-05,6.54964e-08,-2.2169e-11,72653.1,29.6727], Tmin=(100,'K'), Tmax=(879.587,'K')), NASAPolynomial(coeffs=[4.27154,0.0415082,-1.69475e-05,3.00103e-09,-1.98411e-13,72555.1,17.7452], Tmin=(879.587,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(603.275,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(361.68,'J/(mol*K)'), comment="""Thermo group additivity estimation: group(Cs-CsCsCsH) + group(Cs-CsCsHH) + longDistanceInteraction_noncyclic(CsCs-ST) + group(Cs-CsCsHH) + group(Cs-CsHHH) + group(Cs-CsHHH) + group(Cs-CsHHH) + radical(RCCJ) + radical(Tertalkyl) + radical(Isobutyl) + radical(Isobutyl)"""),
)

species(
    label = '[CH2]C([CH2])[CH][CH]C(269)',
    structure = SMILES('[CH2]C([CH2])[CH][CH]C'),
    E0 = (601.594,'kJ/mol'),
    modes = [
        HarmonicOscillator(frequencies=([3000,3050,390,425,1340,1360,335,370,1380,1390,370,380,2900,435,3000,3033.33,3066.67,3100,415,465,780,850,1435,1475,900,1100,2750,2800,2850,1350,1500,750,1050,1375,1000,1338.39,1338.68],'cm^-1')),
        HinderedRotor(inertia=(0.166114,'amu*angstrom^2'), symmetry=1, barrier=(3.81928,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.00300797,'amu*angstrom^2'), symmetry=1, barrier=(3.82284,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.166106,'amu*angstrom^2'), symmetry=1, barrier=(3.8191,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.00300765,'amu*angstrom^2'), symmetry=1, barrier=(3.82471,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.166039,'amu*angstrom^2'), symmetry=1, barrier=(3.81757,'kJ/mol'), semiclassical=False),
    ],
    spinMultiplicity = 5,
    opticalIsomers = 1,
    molecularWeight = (82.1436,'amu'),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[1.51845,0.0570929,-5.08685e-05,3.46663e-08,-1.09234e-11,72442,30.8851], Tmin=(100,'K'), Tmax=(846.196,'K')), NASAPolynomial(coeffs=[4.26798,0.040843,-1.62975e-05,2.88724e-09,-1.92544e-13,72093.1,18.7675], Tmin=(846.196,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(601.594,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(361.68,'J/(mol*K)'), comment="""Thermo group additivity estimation: group(Cs-CsCsCsH) + group(Cs-CsCsHH) + longDistanceInteraction_noncyclic(CsCs-ST) + group(Cs-CsCsHH) + group(Cs-CsHHH) + group(Cs-CsHHH) + group(Cs-CsHHH) + radical(RCCJC) + radical(Isobutyl) + radical(Isobutyl) + radical(Cs_S)"""),
)

species(
    label = '[CH2][C]([CH2])C[CH]C(268)',
    structure = SMILES('[CH2][C]([CH2])C[CH]C'),
    E0 = (592.475,'kJ/mol'),
    modes = [
        HarmonicOscillator(frequencies=([3000,3033.33,3066.67,3100,415,465,780,850,1435,1475,900,1100,2750,2850,1437.5,1250,1305,750,350,2750,2800,2850,1350,1500,750,1050,1375,1000,3025,407.5,1350,352.5,360,370,350,2032.5,2032.73],'cm^-1')),
        HinderedRotor(inertia=(0.213072,'amu*angstrom^2'), symmetry=1, barrier=(4.89894,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.00167157,'amu*angstrom^2'), symmetry=1, barrier=(4.90183,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.0218585,'amu*angstrom^2'), symmetry=1, barrier=(64.075,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.213361,'amu*angstrom^2'), symmetry=1, barrier=(4.90559,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.213085,'amu*angstrom^2'), symmetry=1, barrier=(4.89925,'kJ/mol'), semiclassical=False),
    ],
    spinMultiplicity = 5,
    opticalIsomers = 1,
    molecularWeight = (82.1436,'amu'),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[1.37822,0.066186,-8.8108e-05,8.09903e-08,-2.9222e-11,71344.5,29.9646], Tmin=(100,'K'), Tmax=(894.748,'K')), NASAPolynomial(coeffs=[1.216,0.0459564,-1.90646e-05,3.37205e-09,-2.21395e-13,72212.4,35.4165], Tmin=(894.748,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(592.475,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(361.68,'J/(mol*K)'), comment="""Thermo group additivity estimation: group(Cs-CsCsCsH) + group(Cs-CsCsHH) + longDistanceInteraction_noncyclic(CsCs-ST) + group(Cs-CsCsHH) + group(Cs-CsHHH) + group(Cs-CsHHH) + group(Cs-CsHHH) + radical(Tertalkyl) + radical(Isobutyl) + radical(Isobutyl) + radical(RCCJC)"""),
)

species(
    label = '[CH2][CH]C[C]([CH2])[CH2](296)',
    structure = SMILES('[CH2][CH]C[C]([CH2])[CH2]'),
    E0 = (797.721,'kJ/mol'),
    modes = [
        HarmonicOscillator(frequencies=([2750,2850,1437.5,1250,1305,750,350,360,370,350,3000,3020,3040,3060,3080,3100,415,440,465,780,815,850,1435,1455,1475,900,1000,1100,3025,407.5,1350,352.5,180,2200.61],'cm^-1')),
        HinderedRotor(inertia=(0.00142474,'amu*angstrom^2'), symmetry=1, barrier=(4.89464,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.0188519,'amu*angstrom^2'), symmetry=1, barrier=(64.6986,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.0188544,'amu*angstrom^2'), symmetry=1, barrier=(64.7303,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.213221,'amu*angstrom^2'), symmetry=1, barrier=(4.90236,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.00142194,'amu*angstrom^2'), symmetry=1, barrier=(4.8848,'kJ/mol'), semiclassical=False),
    ],
    spinMultiplicity = 6,
    opticalIsomers = 1,
    molecularWeight = (81.1357,'amu'),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[1.41035,0.0669819,-9.87862e-05,9.35427e-08,-3.37861e-11,96027.4,31.6356], Tmin=(100,'K'), Tmax=(901.877,'K')), NASAPolynomial(coeffs=[1.02444,0.0436355,-1.82801e-05,3.22553e-09,-2.10454e-13,97116.1,39.1074], Tmin=(901.877,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(797.721,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(336.736,'J/(mol*K)'), comment="""Thermo group additivity estimation: group(Cs-CsCsCsH) + group(Cs-CsCsHH) + longDistanceInteraction_noncyclic(CsCs-ST) + group(Cs-CsCsHH) + group(Cs-CsHHH) + group(Cs-CsHHH) + group(Cs-CsHHH) + radical(RCCJ) + radical(Isobutyl) + radical(Tertalkyl) + radical(RCCJC) + radical(Isobutyl)"""),
)

species(
    label = '[CH2][CH][CH]C([CH2])[CH2](297)',
    structure = SMILES('[CH2][CH][CH]C([CH2])[CH2]'),
    E0 = (806.841,'kJ/mol'),
    modes = [
        HarmonicOscillator(frequencies=([3000,3050,390,425,1340,1360,335,370,1380,1390,370,380,2900,435,3000,3020,3040,3060,3080,3100,415,440,465,780,815,850,1435,1455,1475,900,1000,1100,180,1304.94],'cm^-1')),
        HinderedRotor(inertia=(0.00318954,'amu*angstrom^2'), symmetry=1, barrier=(3.85281,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.00318933,'amu*angstrom^2'), symmetry=1, barrier=(3.85326,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.0031898,'amu*angstrom^2'), symmetry=1, barrier=(3.85363,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.167733,'amu*angstrom^2'), symmetry=1, barrier=(3.85651,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.0031882,'amu*angstrom^2'), symmetry=1, barrier=(3.853,'kJ/mol'), semiclassical=False),
    ],
    spinMultiplicity = 6,
    opticalIsomers = 1,
    molecularWeight = (81.1357,'amu'),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[1.52594,0.0581989,-6.27225e-05,4.88731e-08,-1.62498e-11,97126,32.6432], Tmin=(100,'K'), Tmax=(876.166,'K')), NASAPolynomial(coeffs=[4.15784,0.0383785,-1.5428e-05,2.72026e-09,-1.79881e-13,96964.3,22.0033], Tmin=(876.166,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(806.841,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(336.736,'J/(mol*K)'), comment="""Thermo group additivity estimation: group(Cs-CsCsCsH) + group(Cs-CsCsHH) + longDistanceInteraction_noncyclic(CsCs-ST) + group(Cs-CsCsHH) + group(Cs-CsHHH) + group(Cs-CsHHH) + group(Cs-CsHHH) + radical(Cs_S) + radical(Isobutyl) + radical(RCCJ) + radical(Isobutyl) + radical(RCCJC)"""),
)

species(
    label = '[CH2][CH]CC(=C)C(250)',
    structure = SMILES('[CH2][CH]CC(=C)C'),
    E0 = (320.336,'kJ/mol'),
    modes = [
        HarmonicOscillator(frequencies=([3025,407.5,1350,352.5,350,440,435,1725,2950,3100,1380,975,1025,1650,2750,2800,2850,1350,1500,750,1050,1375,1000,2750,2850,1437.5,1250,1305,750,350,3000,3100,440,815,1455,1000,180,1599.32],'cm^-1')),
        HinderedRotor(inertia=(0.0890699,'amu*angstrom^2'), symmetry=1, barrier=(2.04789,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.0890716,'amu*angstrom^2'), symmetry=1, barrier=(2.04793,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.0890702,'amu*angstrom^2'), symmetry=1, barrier=(2.0479,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.0890707,'amu*angstrom^2'), symmetry=1, barrier=(2.04791,'kJ/mol'), semiclassical=False),
    ],
    spinMultiplicity = 3,
    opticalIsomers = 1,
    molecularWeight = (82.1436,'amu'),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[1.82264,0.0512649,-2.77872e-05,7.49961e-09,-8.56105e-13,38601.8,25.3924], Tmin=(100,'K'), Tmax=(1779.92,'K')), NASAPolynomial(coeffs=[8.62855,0.0359701,-1.48978e-05,2.67192e-09,-1.7803e-13,36179,-11.3662], Tmin=(1779.92,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(320.336,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(365.837,'J/(mol*K)'), comment="""Thermo group additivity estimation: group(Cs-CsCsHH) + group(Cs-(Cds-Cds)CsHH) + group(Cs-CsHHH) + group(Cs-(Cds-Cds)HHH) + group(Cds-CdsCsCs) + group(Cds-CdsHH) + radical(RCCJC) + radical(RCCJ)"""),
)

species(
    label = '[CH2]CCC([CH2])=C(298)',
    structure = SMILES('[CH2]CCC([CH2])=C'),
    E0 = (277.388,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
    molecularWeight = (82.1436,'amu'),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[1.08441,0.0545992,-1.97898e-05,-9.49084e-09,6.52881e-12,33475.2,24.9205], Tmin=(100,'K'), Tmax=(1058.22,'K')), NASAPolynomial(coeffs=[12.3682,0.0305284,-1.20085e-05,2.19977e-09,-1.53008e-13,30046.7,-35.0715], Tmin=(1058.22,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(277.388,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(365.837,'J/(mol*K)'), comment="""Thermo group additivity estimation: group(Cs-CsCsHH) + group(Cs-(Cds-Cds)CsHH) + group(Cs-CsHHH) + group(Cs-(Cds-Cds)HHH) + group(Cds-CdsCsCs) + group(Cds-CdsHH) + radical(Allyl_P) + radical(RCCJ)"""),
)

species(
    label = '[CH2]C=CC([CH2])C(251)',
    structure = SMILES('[CH2]C=CC([CH2])C'),
    E0 = (272.282,'kJ/mol'),
    modes = [
        HarmonicOscillator(frequencies=([2995,3025,975,1000,1300,1375,400,500,1630,1680,1380,1390,370,380,2900,435,3000,3033.33,3066.67,3100,415,465,780,850,1435,1475,900,1100,2750,2800,2850,1350,1500,750,1050,1375,1000,180],'cm^-1')),
        HinderedRotor(inertia=(0.0050968,'amu*angstrom^2'), symmetry=1, barrier=(4.60502,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.680192,'amu*angstrom^2'), symmetry=1, barrier=(15.639,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.0173784,'amu*angstrom^2'), symmetry=1, barrier=(15.6342,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(3.30743,'amu*angstrom^2'), symmetry=1, barrier=(76.0443,'kJ/mol'), semiclassical=False),
    ],
    spinMultiplicity = 3,
    opticalIsomers = 1,
    molecularWeight = (82.1436,'amu'),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[1.32272,0.0477992,1.34892e-06,-3.42469e-08,1.64572e-11,32854.3,25.5548], Tmin=(100,'K'), Tmax=(967.428,'K')), NASAPolynomial(coeffs=[12.2803,0.0288316,-1.00792e-05,1.77004e-09,-1.22479e-13,29501.6,-33.3163], Tmin=(967.428,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(272.282,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(365.837,'J/(mol*K)'), comment="""Thermo group additivity estimation: group(Cs-(Cds-Cds)CsCsH) + group(Cs-CsHHH) + group(Cs-CsHHH) + group(Cs-(Cds-Cds)HHH) + group(Cds-CdsCsH) + group(Cds-CdsCsH) + radical(Isobutyl) + radical(Allyl_P)"""),
)

species(
    label = '[CH2][CH]CC[CH][CH2](213)',
    structure = SMILES('[CH2][CH]CC[CH][CH2]'),
    E0 = (607.521,'kJ/mol'),
    modes = [
        HarmonicOscillator(frequencies=([2750,2783.33,2816.67,2850,1425,1450,1225,1275,1270,1340,700,800,300,400,3000,3033.33,3066.67,3100,415,465,780,850,1435,1475,900,1100,3000,3050,390,425,1340,1360,335,370,406.91,2706.24,4000],'cm^-1')),
        HinderedRotor(inertia=(0.0168246,'amu*angstrom^2'), symmetry=1, barrier=(3.0235,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.0168246,'amu*angstrom^2'), symmetry=1, barrier=(3.0235,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.0168246,'amu*angstrom^2'), symmetry=1, barrier=(3.0235,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.0168246,'amu*angstrom^2'), symmetry=1, barrier=(3.0235,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.0168246,'amu*angstrom^2'), symmetry=1, barrier=(3.0235,'kJ/mol'), semiclassical=False),
    ],
    spinMultiplicity = 5,
    opticalIsomers = 1,
    molecularWeight = (82.1436,'amu'),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[2.96636,0.044188,-1.8776e-05,3.20922e-09,-1.84569e-13,73082.1,26.0156], Tmin=(100,'K'), Tmax=(2870.19,'K')), NASAPolynomial(coeffs=[37.1372,0.00431509,-1.98755e-06,2.5034e-10,-8.77445e-15,50274.9,-180.428], Tmin=(2870.19,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(607.521,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(361.68,'J/(mol*K)'), comment="""Thermo group additivity estimation: group(Cs-CsCsHH) + group(Cs-CsCsHH) + group(Cs-CsCsHH) + group(Cs-CsCsHH) + group(Cs-CsHHH) + group(Cs-CsHHH) + radical(RCCJ) + radical(RCCJC) + radical(RCCJ) + radical(RCCJC)"""),
)

species(
    label = '[CH2]C([CH2])C([CH2])[CH2](215)',
    structure = SMILES('[CH2]C([CH2])C([CH2])[CH2]'),
    E0 = (617.076,'kJ/mol'),
    modes = [
        HarmonicOscillator(frequencies=([1380,1383.33,1386.67,1390,370,373.333,376.667,380,2800,3000,430,440,3000,3014.29,3028.57,3042.86,3057.14,3071.43,3085.71,3100,415,431.667,448.333,465,780,803.333,826.667,850,1435,1448.33,1461.67,1475,900,966.667,1033.33,1100,709.057],'cm^-1')),
        HinderedRotor(inertia=(0.00454599,'amu*angstrom^2'), symmetry=1, barrier=(1.6607,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.00496034,'amu*angstrom^2'), symmetry=1, barrier=(1.76166,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.0821712,'amu*angstrom^2'), symmetry=1, barrier=(1.88928,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.0729915,'amu*angstrom^2'), symmetry=1, barrier=(1.67822,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.0734962,'amu*angstrom^2'), symmetry=1, barrier=(1.68982,'kJ/mol'), semiclassical=False),
    ],
    spinMultiplicity = 5,
    opticalIsomers = 1,
    molecularWeight = (82.1436,'amu'),
    collisionModel = TransportData(shapeIndex=2, epsilon=(3248.64,'J/mol'), sigma=(6.11673,'angstroms'), dipoleMoment=(0,'C*m'), polarizability=(0,'angstroms^3'), rotrelaxcollnum=0, comment="""Epsilon & sigma estimated with Tc=507.43 K, Pc=32.21 bar (from Joback method)"""),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[0.668908,0.0650716,-5.68346e-05,2.94123e-08,-6.08665e-12,74344,30.5295], Tmin=(100,'K'), Tmax=(1324.26,'K')), NASAPolynomial(coeffs=[11.093,0.0274374,-7.24257e-06,9.40727e-10,-4.98558e-14,72122.3,-20.6529], Tmin=(1324.26,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(617.076,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(361.68,'J/(mol*K)'), comment="""Thermo group additivity estimation: group(Cs-CsCsCsH) + longDistanceInteraction_noncyclic(CsCs-TT) + group(Cs-CsCsCsH) + longDistanceInteraction_noncyclic(CsCs-TT) + group(Cs-CsHHH) + group(Cs-CsHHH) + group(Cs-CsHHH) + group(Cs-CsHHH) + radical(Isobutyl) + radical(Isobutyl) + radical(Isobutyl) + radical(Isobutyl)"""),
)

species(
    label = '[CH2][CH]CC1CC1(299)',
    structure = SMILES('[CH2][CH]CC1CC1'),
    E0 = (361.383,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
    molecularWeight = (82.1436,'amu'),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[1.83852,0.0358891,2.43859e-05,-4.97208e-08,1.96237e-11,43552.4,25.6475], Tmin=(100,'K'), Tmax=(1004.93,'K')), NASAPolynomial(coeffs=[9.5704,0.0328881,-1.25923e-05,2.31301e-09,-1.62796e-13,40595.9,-18.6703], Tmin=(1004.93,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(361.383,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(369.994,'J/(mol*K)'), comment="""Thermo group additivity estimation: group(Cs-CsCsCsH) + group(Cs-CsCsHH) + group(Cs-CsCsHH) + group(Cs-CsCsHH) + group(Cs-CsCsHH) + group(Cs-CsHHH) + ring(Cyclopropane) + radical(RCCJ) + radical(RCCJC)"""),
)

species(
    label = '[CH2]C1CC([CH2])C1(216)',
    structure = SMILES('[CH2]C1CC([CH2])C1'),
    E0 = (358.394,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
    molecularWeight = (82.1436,'amu'),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[1.96633,0.0264424,6.69093e-05,-1.05639e-07,4.30368e-11,43194.9,21.8847], Tmin=(100,'K'), Tmax=(931.06,'K')), NASAPolynomial(coeffs=[12.9664,0.0258016,-7.16198e-06,1.17495e-09,-8.34205e-14,39126,-41.2495], Tmin=(931.06,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(358.394,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(374.151,'J/(mol*K)'), comment="""Thermo group additivity estimation: group(Cs-CsCsCsH) + group(Cs-CsCsCsH) + group(Cs-CsCsHH) + group(Cs-CsCsHH) + group(Cs-CsHHH) + group(Cs-CsHHH) + ring(Cyclobutane) + radical(Isobutyl) + radical(Isobutyl)"""),
)

species(
    label = '[CH2]C([CH2])CC=C(210)',
    structure = SMILES('[CH2]C([CH2])CC=C'),
    E0 = (337.041,'kJ/mol'),
    modes = [
        HarmonicOscillator(frequencies=([2750,2850,1437.5,1250,1305,750,350,2950,3100,1380,975,1025,1650,1380,1390,370,380,2900,435,3000,3033.33,3066.67,3100,415,465,780,850,1435,1475,900,1100,3010,987.5,1337.5,450,1655,277.213,1590.56],'cm^-1')),
        HinderedRotor(inertia=(0.00211201,'amu*angstrom^2'), symmetry=1, barrier=(0.119627,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.00211184,'amu*angstrom^2'), symmetry=1, barrier=(0.119627,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.132871,'amu*angstrom^2'), symmetry=1, barrier=(7.78714,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(1.13179,'amu*angstrom^2'), symmetry=1, barrier=(69.7505,'kJ/mol'), semiclassical=False),
    ],
    spinMultiplicity = 3,
    opticalIsomers = 1,
    molecularWeight = (82.1436,'amu'),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[1.29623,0.0520097,-1.6468e-05,-1.38423e-08,9.25212e-12,40640.7,26.8619], Tmin=(100,'K'), Tmax=(947.416,'K')), NASAPolynomial(coeffs=[10.6915,0.0302258,-1.0292e-05,1.73515e-09,-1.15664e-13,38057.9,-22.1931], Tmin=(947.416,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(337.041,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(365.837,'J/(mol*K)'), comment="""Thermo group additivity estimation: group(Cs-CsCsCsH) + group(Cs-(Cds-Cds)CsHH) + group(Cs-CsHHH) + group(Cs-CsHHH) + group(Cds-CdsCsH) + group(Cds-CdsHH) + radical(Isobutyl) + radical(Isobutyl)"""),
)

species(
    label = '[CH2][CH]C[CH][CH2](233)',
    structure = SMILES('[CH2][CH]C[CH][CH2]'),
    E0 = (631.301,'kJ/mol'),
    modes = [
        HarmonicOscillator(frequencies=([2750,2850,1437.5,1250,1305,750,350,3000,3033.33,3066.67,3100,415,465,780,850,1435,1475,900,1100,3000,3050,390,425,1340,1360,335,370,1412.64,1413.59],'cm^-1')),
        HinderedRotor(inertia=(0.00561943,'amu*angstrom^2'), symmetry=1, barrier=(7.94511,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.00562752,'amu*angstrom^2'), symmetry=1, barrier=(7.94764,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.00561538,'amu*angstrom^2'), symmetry=1, barrier=(7.94369,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.00558985,'amu*angstrom^2'), symmetry=1, barrier=(7.94911,'kJ/mol'), semiclassical=False),
    ],
    spinMultiplicity = 5,
    opticalIsomers = 1,
    molecularWeight = (68.117,'amu'),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[2.26905,0.0435843,-4.87511e-05,4.66777e-08,-1.86946e-11,75984.9,26.4118], Tmin=(100,'K'), Tmax=(836.611,'K')), NASAPolynomial(coeffs=[0.578571,0.0388266,-1.7199e-05,3.1897e-09,-2.17338e-13,76717.1,36.9514], Tmin=(836.611,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(631.301,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(291.007,'J/(mol*K)'), comment="""Thermo group additivity estimation: group(Cs-CsCsHH) + group(Cs-CsCsHH) + group(Cs-CsCsHH) + group(Cs-CsHHH) + group(Cs-CsHHH) + radical(RCCJC) + radical(RCCJ) + radical(RCCJ) + radical(RCCJC)"""),
)

species(
    label = '[CH][CH2](235)',
    structure = SMILES('[CH][CH2]'),
    E0 = (556.742,'kJ/mol'),
    modes = [
        HarmonicOscillator(frequencies=([3000,3100,440,815,1455,1000,1101.59,1101.66],'cm^-1')),
        HinderedRotor(inertia=(0.00420681,'amu*angstrom^2'), symmetry=1, barrier=(3.6236,'kJ/mol'), semiclassical=False),
    ],
    spinMultiplicity = 4,
    opticalIsomers = 1,
    molecularWeight = (27.0452,'amu'),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[3.77499,-0.000463351,3.18197e-05,-4.30828e-08,1.77627e-11,66973.8,8.7898], Tmin=(100,'K'), Tmax=(870.332,'K')), NASAPolynomial(coeffs=[6.06982,0.00332463,5.85313e-07,-2.32963e-10,1.82425e-14,66031.4,-5.08173], Tmin=(870.332,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(556.742,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(103.931,'J/(mol*K)'), comment="""Thermo group additivity estimation: group(Cs-CsHHH) + group(Cs-CsHHH) + radical(CCJ) + radical(CCJ2_triplet)"""),
)

species(
    label = '[CH2]C([CH2])[CH2](226)',
    structure = SMILES('[CH2]C([CH2])[CH2]'),
    E0 = (461.901,'kJ/mol'),
    modes = [
        HarmonicOscillator(frequencies=([1380,1390,370,380,2900,435,3000,3020,3040,3060,3080,3100,415,440,465,780,815,850,1435,1455,1475,900,1000,1100],'cm^-1')),
        HinderedRotor(inertia=(0.00529301,'amu*angstrom^2'), symmetry=1, barrier=(3.08377,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.00523118,'amu*angstrom^2'), symmetry=1, barrier=(3.03916,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.130467,'amu*angstrom^2'), symmetry=1, barrier=(2.9997,'kJ/mol'), semiclassical=False),
    ],
    spinMultiplicity = 4,
    opticalIsomers = 1,
    molecularWeight = (55.0984,'amu'),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[1.9185,0.0363463,-2.61059e-05,1.15462e-08,-2.04218e-12,55637,20.3501], Tmin=(100,'K'), Tmax=(1631.3,'K')), NASAPolynomial(coeffs=[7.087,0.0180707,-4.14987e-06,4.68151e-10,-2.1815e-14,54696.1,-4.82954], Tmin=(1631.3,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(461.901,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(245.277,'J/(mol*K)'), comment="""Thermo group additivity estimation: group(Cs-CsCsCsH) + group(Cs-CsHHH) + group(Cs-CsHHH) + group(Cs-CsHHH) + radical(Isobutyl) + radical(Isobutyl) + radical(Isobutyl)"""),
)

species(
    label = '[CH2][C]CC([CH2])[CH2](300)',
    structure = SMILES('[CH2][C]CC([CH2])[CH2]'),
    E0 = (866.067,'kJ/mol'),
    modes = [
        HarmonicOscillator(frequencies=([2750,2850,1437.5,1250,1305,750,350,1380,1390,370,380,2900,435,3000,3020,3040,3060,3080,3100,415,440,465,780,815,850,1435,1455,1475,900,1000,1100,200,800,1600],'cm^-1')),
        HinderedRotor(inertia=(0.156089,'amu*angstrom^2'), symmetry=1, barrier=(3.5888,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.156089,'amu*angstrom^2'), symmetry=1, barrier=(3.5888,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.156089,'amu*angstrom^2'), symmetry=1, barrier=(3.5888,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.156089,'amu*angstrom^2'), symmetry=1, barrier=(3.5888,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.156089,'amu*angstrom^2'), symmetry=1, barrier=(3.5888,'kJ/mol'), semiclassical=False),
    ],
    spinMultiplicity = 6,
    opticalIsomers = 1,
    molecularWeight = (81.1357,'amu'),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[0.928002,0.064135,-5.80072e-05,3.03055e-08,-6.51296e-12,104277,29.4001], Tmin=(100,'K'), Tmax=(1118.32,'K')), NASAPolynomial(coeffs=[10.949,0.0282919,-9.9306e-06,1.64537e-09,-1.05952e-13,102036,-20.0658], Tmin=(1118.32,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(866.067,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(336.736,'J/(mol*K)'), comment="""Thermo group additivity estimation: group(Cs-CsCsCsH) + group(Cs-CsCsHH) + longDistanceInteraction_noncyclic(CsCs-ST) + group(Cs-CsCsHH) + group(Cs-CsHHH) + group(Cs-CsHHH) + group(Cs-CsHHH) + radical(RCCJ) + radical(Isobutyl) + radical(CCJ2_triplet) + radical(Isobutyl)"""),
)

species(
    label = '[CH]C([CH2])C[CH][CH2](301)',
    structure = SMILES('[CH]C([CH2])C[CH][CH2]'),
    E0 = (855.431,'kJ/mol'),
    modes = [
        HarmonicOscillator(frequencies=([2750,2850,1437.5,1250,1305,750,350,1380,1390,370,380,2900,435,3000,3033.33,3066.67,3100,415,465,780,850,1435,1475,900,1100,3025,407.5,1350,352.5,210.214,429.066,1087.85,1965.66,3037.96],'cm^-1')),
        HinderedRotor(inertia=(0.0721528,'amu*angstrom^2'), symmetry=1, barrier=(1.66656,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.0721528,'amu*angstrom^2'), symmetry=1, barrier=(1.66656,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.0721528,'amu*angstrom^2'), symmetry=1, barrier=(1.66656,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.0721528,'amu*angstrom^2'), symmetry=1, barrier=(1.66656,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.0721528,'amu*angstrom^2'), symmetry=1, barrier=(1.66656,'kJ/mol'), semiclassical=False),
    ],
    spinMultiplicity = 6,
    opticalIsomers = 1,
    molecularWeight = (81.1357,'amu'),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[1.3854,0.0592901,-5.27553e-05,3.01056e-08,-7.57637e-12,102977,30.5122], Tmin=(100,'K'), Tmax=(929.4,'K')), NASAPolynomial(coeffs=[6.90797,0.0355223,-1.43962e-05,2.59096e-09,-1.75335e-13,101951,4.27322], Tmin=(929.4,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(855.431,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(336.736,'J/(mol*K)'), comment="""Thermo group additivity estimation: group(Cs-CsCsCsH) + group(Cs-CsCsHH) + longDistanceInteraction_noncyclic(CsCs-ST) + group(Cs-CsCsHH) + group(Cs-CsHHH) + group(Cs-CsHHH) + group(Cs-CsHHH) + radical(RCCJC) + radical(RCCJ) + radical(CCJ2_triplet) + radical(Isobutyl)"""),
)

species(
    label = '[CH][CH]CC([CH2])[CH2](302)',
    structure = SMILES('[CH][CH]CC([CH2])[CH2]'),
    E0 = (855.267,'kJ/mol'),
    modes = [
        HarmonicOscillator(frequencies=([2750,2850,1437.5,1250,1305,750,350,1380,1390,370,380,2900,435,3000,3033.33,3066.67,3100,415,465,780,850,1435,1475,900,1100,3025,407.5,1350,352.5,220.196,1047.78,1397.03,1690.21,2212.38],'cm^-1')),
        HinderedRotor(inertia=(0.11319,'amu*angstrom^2'), symmetry=1, barrier=(2.78752,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.11319,'amu*angstrom^2'), symmetry=1, barrier=(2.78752,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.11319,'amu*angstrom^2'), symmetry=1, barrier=(2.78752,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.11319,'amu*angstrom^2'), symmetry=1, barrier=(2.78752,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.11319,'amu*angstrom^2'), symmetry=1, barrier=(2.78752,'kJ/mol'), semiclassical=False),
    ],
    spinMultiplicity = 6,
    opticalIsomers = 1,
    molecularWeight = (81.1357,'amu'),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[1.52843,0.0577316,-4.24667e-05,7.26199e-09,7.53331e-12,102951,29.3345], Tmin=(100,'K'), Tmax=(654.245,'K')), NASAPolynomial(coeffs=[7.16544,0.0340336,-1.28177e-05,2.20259e-09,-1.44655e-13,101983,2.76987], Tmin=(654.245,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(855.267,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(336.736,'J/(mol*K)'), comment="""Thermo group additivity estimation: group(Cs-CsCsCsH) + group(Cs-CsCsHH) + longDistanceInteraction_noncyclic(CsCs-ST) + group(Cs-CsCsHH) + group(Cs-CsHHH) + group(Cs-CsHHH) + group(Cs-CsHHH) + radical(RCCJC) + radical(CCJ2_triplet) + radical(Isobutyl) + radical(Isobutyl)"""),
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
    E0 = (612.298,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS2',
    E0 = (694.928,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS3',
    E0 = (696.943,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS4',
    E0 = (740.761,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS5',
    E0 = (662.428,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS6',
    E0 = (753.786,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS7',
    E0 = (749.398,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS8',
    E0 = (730.238,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS9',
    E0 = (728.434,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS10',
    E0 = (756.629,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS11',
    E0 = (743.424,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS12',
    E0 = (968.585,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS13',
    E0 = (1009.53,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS14',
    E0 = (1018.65,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS15',
    E0 = (635.16,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS16',
    E0 = (675.699,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS17',
    E0 = (675.699,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS18',
    E0 = (769.617,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS19',
    E0 = (774.395,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS20',
    E0 = (620.388,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS21',
    E0 = (620.583,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS22',
    E0 = (612.298,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS23',
    E0 = (1046.99,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS24',
    E0 = (1052.96,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS25',
    E0 = (1077.87,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS26',
    E0 = (1067.24,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS27',
    E0 = (1067.07,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

reaction(
    label = 'reaction1',
    reactants = ['[CH2][CH]CC([CH2])[CH2](214)'],
    products = ['[CH2]C=C(66)', '[CH2]C=C(66)'],
    transitionState = 'TS1',
    kinetics = Arrhenius(A=(5e+12,'s^-1'), n=0, Ea=(0,'kJ/mol'), T0=(1,'K'), Tmin=(300,'K'), Tmax=(1500,'K'), comment="""Exact match found for rate rule [RJJ]
Euclidian distance = 0
family: 1,4_Linear_birad_scission"""),
)

reaction(
    label = 'reaction2',
    reactants = ['H(7)', '[CH2][CH]CC([CH2])=C(292)'],
    products = ['[CH2][CH]CC([CH2])[CH2](214)'],
    transitionState = 'TS2',
    kinetics = Arrhenius(A=(170.395,'m^3/(mol*s)'), n=1.5621, Ea=(11.2886,'kJ/mol'), T0=(1,'K'), comment="""Estimated using an average for rate rule [Cds_Cds;HJ]
Euclidian distance = 0
family: R_Addition_MultipleBond"""),
)

reaction(
    label = 'reaction3',
    reactants = ['H(7)', '[CH2]C=CC([CH2])[CH2](293)'],
    products = ['[CH2][CH]CC([CH2])[CH2](214)'],
    transitionState = 'TS3',
    kinetics = Arrhenius(A=(8.22e+08,'cm^3/(mol*s)'), n=1.533, Ea=(7.77387,'kJ/mol'), T0=(1,'K'), comment="""From training reaction 192 used for Cd_R;HJ
Exact match found for rate rule [Cd_R;HJ]
Euclidian distance = 0
family: R_Addition_MultipleBond"""),
)

reaction(
    label = 'reaction4',
    reactants = ['[CH2][CH]CC=C(188)', 'CH2(T)(26)'],
    products = ['[CH2][CH]CC([CH2])[CH2](214)'],
    transitionState = 'TS4',
    kinetics = Arrhenius(A=(4.17e+07,'cm^3/(mol*s)'), n=1.64, Ea=(0,'kJ/mol'), T0=(1,'K'), Tmin=(300,'K'), Tmax=(1500,'K'), comment="""Estimated using template [Cds-CsH_Cds-HH;Y_1centerbirad] for rate rule [Cds-CsH_Cds-HH;CH2_triplet]
Euclidian distance = 1.0
family: R_Addition_MultipleBond
Ea raised from -5.9 to 0 kJ/mol."""),
)

reaction(
    label = 'reaction5',
    reactants = ['[CH2][CH][CH2](219)', '[CH2]C=C(66)'],
    products = ['[CH2][CH]CC([CH2])[CH2](214)'],
    transitionState = 'TS5',
    kinetics = Arrhenius(A=(0.00359021,'m^3/(mol*s)'), n=2.50446, Ea=(20.5128,'kJ/mol'), T0=(1,'K'), comment="""Estimated using an average for rate rule [Cds_Cds;CJ]
Euclidian distance = 0
Multiplied by reaction path degeneracy 2.0
family: R_Addition_MultipleBond"""),
)

reaction(
    label = 'reaction6',
    reactants = ['[CH2][CH]CC([CH2])[CH2](214)'],
    products = ['[CH2][CH]C[C]([CH2])C(257)'],
    transitionState = 'TS6',
    kinetics = Arrhenius(A=(34962.3,'s^-1'), n=2.56136, Ea=(141.488,'kJ/mol'), T0=(1,'K'), comment="""Estimated using an average for rate rule [R2H_S;C_rad_out_2H;XH_out]
Euclidian distance = 0
Multiplied by reaction path degeneracy 2.0
family: intra_H_migration"""),
)

reaction(
    label = 'reaction7',
    reactants = ['[CH2]C[CH]C([CH2])[CH2](294)'],
    products = ['[CH2][CH]CC([CH2])[CH2](214)'],
    transitionState = 'TS7',
    kinetics = Arrhenius(A=(1.91383e+06,'s^-1'), n=2.02407, Ea=(137.004,'kJ/mol'), T0=(1,'K'), comment="""Estimated using an average for rate rule [R2H_S;C_rad_out_H/NonDeC;XH_out]
Euclidian distance = 0
Multiplied by reaction path degeneracy 2.0
family: intra_H_migration"""),
)

reaction(
    label = 'reaction8',
    reactants = ['[CH2][CH]CC([CH2])[CH2](214)'],
    products = ['[CH2][CH][CH]C([CH2])C(258)'],
    transitionState = 'TS8',
    kinetics = Arrhenius(A=(333380,'s^-1'), n=2.17519, Ea=(117.939,'kJ/mol'), T0=(1,'K'), comment="""Estimated using an average for rate rule [R3H_SS_Cs;C_rad_out_2H;XH_out]
Euclidian distance = 0
Multiplied by reaction path degeneracy 4.0
family: intra_H_migration"""),
)

reaction(
    label = 'reaction9',
    reactants = ['[CH2][CH]CC([CH2])[CH2](214)'],
    products = ['[CH2]CC[C]([CH2])[CH2](295)'],
    transitionState = 'TS9',
    kinetics = Arrhenius(A=(285601,'s^-1'), n=2.01653, Ea=(116.136,'kJ/mol'), T0=(1,'K'), comment="""Estimated using an average for rate rule [R3H_SS_Cs;Y_rad_out;XH_out]
Euclidian distance = 0
family: intra_H_migration"""),
)

reaction(
    label = 'reaction10',
    reactants = ['[CH2][CH]CC([CH2])[CH2](214)'],
    products = ['[CH2]C([CH2])[CH][CH]C(269)'],
    transitionState = 'TS10',
    kinetics = Arrhenius(A=(2.50921e+06,'s^-1'), n=1.8375, Ea=(144.331,'kJ/mol'), T0=(1,'K'), comment="""Estimated using template [R3Hall;C_rad_out_2H;Cs_H_out_H/NonDeC] for rate rule [R3HJ;C_rad_out_2H;Cs_H_out_H/NonDeC]
Euclidian distance = 1.0
Multiplied by reaction path degeneracy 2.0
family: intra_H_migration"""),
)

reaction(
    label = 'reaction11',
    reactants = ['[CH2][CH]CC([CH2])[CH2](214)'],
    products = ['[CH2][C]([CH2])C[CH]C(268)'],
    transitionState = 'TS11',
    kinetics = Arrhenius(A=(434.916,'s^-1'), n=2.67444, Ea=(131.125,'kJ/mol'), T0=(1,'K'), comment="""Estimated using template [R4Hall;C_rad_out_2H;XH_out] for rate rule [R4HJ_1;C_rad_out_2H;XH_out]
Euclidian distance = 1.0
family: intra_H_migration"""),
)

reaction(
    label = 'reaction12',
    reactants = ['[CH2][CH][CH2](219)', '[CH2][CH][CH2](219)'],
    products = ['[CH2][CH]CC([CH2])[CH2](214)'],
    transitionState = 'TS12',
    kinetics = Arrhenius(A=(3.9578e+07,'m^3/(mol*s)'), n=-0.126319, Ea=(0,'kJ/mol'), T0=(1,'K'), comment="""Estimated using an average for rate rule [Y_rad;Y_rad]
Euclidian distance = 0
Multiplied by reaction path degeneracy 2.0
family: R_Recombination
Ea raised from -15.6 to -15.6 kJ/mol.
Ea raised from -15.6 to 0 kJ/mol."""),
)

reaction(
    label = 'reaction13',
    reactants = ['[CH2][CH]C[C]([CH2])[CH2](296)', 'H(7)'],
    products = ['[CH2][CH]CC([CH2])[CH2](214)'],
    transitionState = 'TS13',
    kinetics = Arrhenius(A=(4.34078e+06,'m^3/(mol*s)'), n=0.278577, Ea=(0,'kJ/mol'), T0=(1,'K'), comment="""Estimated using an average for rate rule [Y_rad;H_rad]
Euclidian distance = 0
family: R_Recombination
Ea raised from -1.4 to 0 kJ/mol."""),
)

reaction(
    label = 'reaction14',
    reactants = ['[CH2][CH][CH]C([CH2])[CH2](297)', 'H(7)'],
    products = ['[CH2][CH]CC([CH2])[CH2](214)'],
    transitionState = 'TS14',
    kinetics = Arrhenius(A=(4.34078e+06,'m^3/(mol*s)'), n=0.278577, Ea=(0,'kJ/mol'), T0=(1,'K'), comment="""Estimated using an average for rate rule [Y_rad;H_rad]
Euclidian distance = 0
family: R_Recombination
Ea raised from -1.4 to 0 kJ/mol."""),
)

reaction(
    label = 'reaction15',
    reactants = ['[CH2][CH]CC([CH2])[CH2](214)'],
    products = ['[CH2][CH]CC(=C)C(250)'],
    transitionState = 'TS15',
    kinetics = Arrhenius(A=(3.898e+11,'s^-1'), n=0.486, Ea=(22.8614,'kJ/mol'), T0=(1,'K'), comment="""Estimated using an average for rate rule [R2radExo;Y_rad;XH_Rrad]
Euclidian distance = 0
Multiplied by reaction path degeneracy 2.0
family: Intra_Disproportionation"""),
)

reaction(
    label = 'reaction16',
    reactants = ['[CH2][CH]CC([CH2])[CH2](214)'],
    products = ['[CH2]CCC([CH2])=C(298)'],
    transitionState = 'TS16',
    kinetics = Arrhenius(A=(7.437e+08,'s^-1'), n=1.045, Ea=(63.4002,'kJ/mol'), T0=(1,'K'), comment="""Estimated using an average for rate rule [R3radExo;Y_rad;XH_Rrad]
Euclidian distance = 0
family: Intra_Disproportionation"""),
)

reaction(
    label = 'reaction17',
    reactants = ['[CH2][CH]CC([CH2])[CH2](214)'],
    products = ['[CH2]C=CC([CH2])C(251)'],
    transitionState = 'TS17',
    kinetics = Arrhenius(A=(2.9748e+09,'s^-1'), n=1.045, Ea=(63.4002,'kJ/mol'), T0=(1,'K'), comment="""Estimated using an average for rate rule [R3radExo;Y_rad;XH_Rrad]
Euclidian distance = 0
Multiplied by reaction path degeneracy 4.0
family: Intra_Disproportionation"""),
)

reaction(
    label = 'reaction18',
    reactants = ['[CH2][CH]CC([CH2])[CH2](214)'],
    products = ['[CH2][CH]CC[CH][CH2](213)'],
    transitionState = 'TS18',
    kinetics = Arrhenius(A=(2.66e+08,'s^-1'), n=1.36, Ea=(157.318,'kJ/mol'), T0=(1,'K'), comment="""Estimated using an average for rate rule [cCsCJ;CsJ-HH;C]
Euclidian distance = 0
Multiplied by reaction path degeneracy 2.0
family: 1,2_shiftC"""),
)

reaction(
    label = 'reaction11',
    reactants = ['[CH2]C([CH2])C([CH2])[CH2](215)'],
    products = ['[CH2][CH]CC([CH2])[CH2](214)'],
    transitionState = 'TS19',
    kinetics = Arrhenius(A=(5.32e+08,'s^-1'), n=1.36, Ea=(157.318,'kJ/mol'), T0=(1,'K'), comment="""Estimated using an average for rate rule [cCsCJ;CsJ-HH;C]
Euclidian distance = 0
Multiplied by reaction path degeneracy 4.0
family: 1,2_shiftC"""),
)

reaction(
    label = 'reaction20',
    reactants = ['[CH2][CH]CC([CH2])[CH2](214)'],
    products = ['[CH2][CH]CC1CC1(299)'],
    transitionState = 'TS20',
    kinetics = Arrhenius(A=(7.38971e+10,'s^-1'), n=0.0476667, Ea=(8.08907,'kJ/mol'), T0=(1,'K'), comment="""Estimated using template [Rn;C_rad_out_2H;Cpri_rad_out_2H] for rate rule [R3_SS;C_rad_out_2H;Cpri_rad_out_2H]
Euclidian distance = 2.0
family: Birad_recombination"""),
)

reaction(
    label = 'reaction21',
    reactants = ['[CH2][CH]CC([CH2])[CH2](214)'],
    products = ['[CH2]C1CC([CH2])C1(216)'],
    transitionState = 'TS21',
    kinetics = Arrhenius(A=(3.24e+12,'s^-1'), n=-0.305, Ea=(8.28432,'kJ/mol'), T0=(1,'K'), Tmin=(600,'K'), Tmax=(2000,'K'), comment="""Estimated using an average for rate rule [R4_SSS;C_rad_out_2H;Ypri_rad_out]
Euclidian distance = 0
Multiplied by reaction path degeneracy 2.0
family: Birad_recombination"""),
)

reaction(
    label = 'reaction22',
    reactants = ['[CH2][CH]CC([CH2])[CH2](214)'],
    products = ['[CH2]C([CH2])CC=C(210)'],
    transitionState = 'TS22',
    kinetics = Arrhenius(A=(6.31e+07,'s^-1'), n=0, Ea=(0,'kJ/mol'), T0=(1,'K'), Tmin=(300,'K'), Tmax=(1500,'K'), comment="""From training reaction 2 used for Y_12_10
Exact match found for rate rule [Y_12_10]
Euclidian distance = 0
family: 1,2-Birad_to_alkene"""),
)

reaction(
    label = 'reaction23',
    reactants = ['[CH2][CH]C[CH][CH2](233)', 'CH2(T)(26)'],
    products = ['[CH2][CH]CC([CH2])[CH2](214)'],
    transitionState = 'TS23',
    kinetics = Arrhenius(A=(2.29708e+06,'m^3/(mol*s)'), n=0.575199, Ea=(34.3157,'kJ/mol'), T0=(1,'K'), comment="""Estimated using an average for rate rule [Y_rad;Birad]
Euclidian distance = 0
Multiplied by reaction path degeneracy 2.0
family: Birad_R_Recombination"""),
)

reaction(
    label = 'reaction24',
    reactants = ['[CH][CH2](235)', '[CH2]C([CH2])[CH2](226)'],
    products = ['[CH2][CH]CC([CH2])[CH2](214)'],
    transitionState = 'TS24',
    kinetics = Arrhenius(A=(3.44562e+06,'m^3/(mol*s)'), n=0.575199, Ea=(34.3157,'kJ/mol'), T0=(1,'K'), comment="""Estimated using template [Y_rad;Birad] for rate rule [C_rad/H2/Cs;Birad]
Euclidian distance = 3.0
Multiplied by reaction path degeneracy 3.0
family: Birad_R_Recombination"""),
)

reaction(
    label = 'reaction25',
    reactants = ['[CH2][C]CC([CH2])[CH2](300)', 'H(7)'],
    products = ['[CH2][CH]CC([CH2])[CH2](214)'],
    transitionState = 'TS25',
    kinetics = Arrhenius(A=(1e+07,'m^3/(mol*s)'), n=0, Ea=(0,'kJ/mol'), T0=(1,'K'), comment="""Estimated using an average for rate rule [H_rad;Birad]
Euclidian distance = 0
family: Birad_R_Recombination"""),
)

reaction(
    label = 'reaction26',
    reactants = ['[CH]C([CH2])C[CH][CH2](301)', 'H(7)'],
    products = ['[CH2][CH]CC([CH2])[CH2](214)'],
    transitionState = 'TS26',
    kinetics = Arrhenius(A=(1e+07,'m^3/(mol*s)'), n=0, Ea=(0,'kJ/mol'), T0=(1,'K'), comment="""Estimated using an average for rate rule [H_rad;Birad]
Euclidian distance = 0
family: Birad_R_Recombination"""),
)

reaction(
    label = 'reaction27',
    reactants = ['[CH][CH]CC([CH2])[CH2](302)', 'H(7)'],
    products = ['[CH2][CH]CC([CH2])[CH2](214)'],
    transitionState = 'TS27',
    kinetics = Arrhenius(A=(1e+07,'m^3/(mol*s)'), n=0, Ea=(0,'kJ/mol'), T0=(1,'K'), comment="""Estimated using an average for rate rule [H_rad;Birad]
Euclidian distance = 0
family: Birad_R_Recombination"""),
)

network(
    label = '122',
    isomers = [
        '[CH2][CH]CC([CH2])[CH2](214)',
    ],
    reactants = [
        ('[CH2]C=C(66)', '[CH2]C=C(66)'),
    ],
    bathGas = {
        'Ne': 0.25,
        'He': 0.25,
        'N2': 0.25,
        'Ar': 0.25,
    },
)

pressureDependence(
    label = '122',
    Tmin = (300,'K'),
    Tmax = (3000,'K'),
    Tcount = 8,
    Tlist = ([302.617,324.619,374.997,470.374,649.057,1000.02,1706.11,2761.25],'K'),
    Pmin = (0.001,'bar'),
    Pmax = (100,'bar'),
    Pcount = 5,
    Plist = ([0.00132544,0.0107284,0.316228,9.32101,75.4469],'bar'),
    maximumGrainSize = (0.5,'kcal/mol'),
    minimumGrainCount = 250,
    method = 'modified strong collision',
    interpolationModel = ('Chebyshev', 6, 4),
    activeKRotor = True,
    activeJRotor = True,
    rmgmode = True,
)

