species(
    label = '[CH]C(=[CH])C([CH2])[CH2](19200)',
    structure = SMILES('[CH]C(=[CH])C([CH2])[CH2]'),
    E0 = (944.8,'kJ/mol'),
    modes = [
        HarmonicOscillator(frequencies=([3000,3033.33,3066.67,3100,415,465,780,850,1435,1475,900,1100,3120,650,792.5,1650,1380,1390,370,380,2900,435,350,440,435,1725,492.573,492.856,493.377],'cm^-1')),
        HinderedRotor(inertia=(0.00069575,'amu*angstrom^2'), symmetry=1, barrier=(0.119627,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.310712,'amu*angstrom^2'), symmetry=1, barrier=(53.542,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.31084,'amu*angstrom^2'), symmetry=1, barrier=(53.541,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.310118,'amu*angstrom^2'), symmetry=1, barrier=(53.5398,'kJ/mol'), semiclassical=False),
    ],
    spinMultiplicity = 6,
    opticalIsomers = 1,
    molecularWeight = (79.1198,'amu'),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[1.16628,0.0586672,-4.15569e-05,1.0902e-08,1.3757e-12,113739,26.6969], Tmin=(100,'K'), Tmax=(889.962,'K')), NASAPolynomial(coeffs=[10.5201,0.0281773,-9.63713e-06,1.57586e-09,-1.01524e-13,111616,-19.9098], Tmin=(889.962,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(944.8,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(291.007,'J/(mol*K)'), comment="""Thermo group additivity estimation: group(Cs-(Cds-Cds)CsCsH) + group(Cs-CsHHH) + group(Cs-CsHHH) + group(Cs-(Cds-Cds)HHH) + group(Cds-CdsCsCs) + longDistanceInteraction_noncyclic(CdCs-ST) + group(Cds-CdsHH) + radical(AllylJ2_triplet) + radical(Isobutyl) + radical(Isobutyl) + radical(Cds_P)"""),
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
    label = '[CH]=C=[CH](18734)',
    structure = SMILES('[CH]=C=[CH]'),
    E0 = (491.681,'kJ/mol'),
    modes = [
        HarmonicOscillator(frequencies=([540,610,2055,239.877,511.233,1743.98,1746.51,1747.6,1753.44],'cm^-1')),
    ],
    spinMultiplicity = 3,
    opticalIsomers = 1,
    molecularWeight = (38.048,'amu'),
    collisionModel = TransportData(shapeIndex=2, epsilon=(1737.73,'J/mol'), sigma=(4.1,'angstroms'), dipoleMoment=(0,'C*m'), polarizability=(0,'angstroms^3'), rotrelaxcollnum=1.0, comment="""GRI-Mech"""),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[2.766,0.0170203,-1.57568e-05,7.95984e-09,-1.4265e-12,59188.9,11.2142], Tmin=(100,'K'), Tmax=(1806.04,'K')), NASAPolynomial(coeffs=[4.81405,0.00509933,2.77647e-07,-2.23082e-10,1.96202e-14,59653.5,3.45727], Tmin=(1806.04,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(491.681,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(108.088,'J/(mol*K)'), comment="""Thermo group additivity estimation: group(Cds-CdsHH) + group(Cds-CdsHH) + group(Cdd-CdsCds) + radical(C=C=CJ) + radical(C=C=CJ)"""),
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
    label = '[CH]C(=[CH])C([CH2])=C(19688)',
    structure = SMILES('[CH]C(=[CH])C([CH2])=C'),
    E0 = (788.02,'kJ/mol'),
    modes = [
        HarmonicOscillator(frequencies=([3120,650,792.5,1650,2950,3100,1380,975,1025,1650,325,375,415,465,420,450,1700,1750,3000,3100,440,815,1455,1000,180,180,180],'cm^-1')),
        HinderedRotor(inertia=(2.12583,'amu*angstrom^2'), symmetry=1, barrier=(48.8771,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(2.12507,'amu*angstrom^2'), symmetry=1, barrier=(48.8596,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(2.12471,'amu*angstrom^2'), symmetry=1, barrier=(48.8512,'kJ/mol'), semiclassical=False),
    ],
    spinMultiplicity = 5,
    opticalIsomers = 1,
    molecularWeight = (78.1118,'amu'),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[0.946564,0.0575781,-3.03186e-05,-6.64489e-09,8.19087e-12,94895.5,20.5622], Tmin=(100,'K'), Tmax=(961.756,'K')), NASAPolynomial(coeffs=[15.1758,0.0205555,-7.13489e-06,1.23977e-09,-8.55912e-14,91133.7,-52.8585], Tmin=(961.756,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(788.02,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(270.22,'J/(mol*K)'), comment="""Thermo group additivity estimation: group(Cs-(Cds-Cds)HHH) + group(Cs-(Cds-Cds)HHH) + group(Cds-Cds(Cds-Cds)Cs) + group(Cds-Cds(Cds-Cds)Cs) + group(Cds-CdsHH) + group(Cds-CdsHH) + radical(AllylJ2_triplet) + radical(Cds_P) + radical(Allyl_P)"""),
)

species(
    label = '[CH](2815)',
    structure = SMILES('[CH]'),
    E0 = (585.033,'kJ/mol'),
    modes = [
        HarmonicOscillator(frequencies=([4000],'cm^-1')),
    ],
    spinMultiplicity = 4,
    opticalIsomers = 1,
    molecularWeight = (13.0186,'amu'),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[4.1763,-0.00339736,5.29655e-06,-3.21799e-09,7.28313e-13,70356.4,-0.99239], Tmin=(100,'K'), Tmax=(1260.74,'K')), NASAPolynomial(coeffs=[3.26554,0.000229807,1.03509e-07,-7.93772e-12,-2.40435e-16,70527.4,3.38009], Tmin=(1260.74,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(585.033,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(108.088,'J/(mol*K)'), comment="""Thermo library: primaryThermoLibrary + radical(CJ3)"""),
)

species(
    label = 'C#CC([CH2])[CH2](6419)',
    structure = SMILES('C#CC([CH2])[CH2]'),
    E0 = (526.841,'kJ/mol'),
    modes = [
        HarmonicOscillator(frequencies=([750,770,3400,2100,1380,1390,370,380,2900,435,2175,525,3000,3033.33,3066.67,3100,415,465,780,850,1435,1475,900,1100],'cm^-1')),
        HinderedRotor(inertia=(3.29698,'amu*angstrom^2'), symmetry=1, barrier=(75.804,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.218602,'amu*angstrom^2'), symmetry=1, barrier=(75.7941,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.0119489,'amu*angstrom^2'), symmetry=1, barrier=(75.837,'kJ/mol'), semiclassical=False),
    ],
    spinMultiplicity = 3,
    opticalIsomers = 1,
    molecularWeight = (66.1011,'amu'),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[1.59865,0.0445704,-4.03451e-05,2.09026e-08,-4.1805e-12,63457.8,21.3678], Tmin=(100,'K'), Tmax=(1443.27,'K')), NASAPolynomial(coeffs=[9.29133,0.0160372,-3.19385e-06,2.79115e-10,-8.33799e-15,61988.6,-15.9649], Tmin=(1443.27,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(526.841,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(245.277,'J/(mol*K)'), comment="""Thermo group additivity estimation: group(Cs-CtCsCsH) + group(Cs-CsHHH) + group(Cs-CsHHH) + group(Ct-CtCs) + group(Ct-CtH) + radical(Isobutyl) + radical(Isobutyl)"""),
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
    label = '[CH]C(=[CH])C=C(19261)',
    structure = SMILES('[CH]C(=[CH])C=C'),
    E0 = (674.111,'kJ/mol'),
    modes = [
        HarmonicOscillator(frequencies=([3120,650,792.5,1650,2950,3100,1380,975,1025,1650,350,440,435,1725,3010,987.5,1337.5,450,1655,180,180,180],'cm^-1')),
        HinderedRotor(inertia=(2.10119,'amu*angstrom^2'), symmetry=1, barrier=(48.3106,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(2.0992,'amu*angstrom^2'), symmetry=1, barrier=(48.2648,'kJ/mol'), semiclassical=False),
    ],
    spinMultiplicity = 4,
    opticalIsomers = 1,
    molecularWeight = (65.0932,'amu'),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[1.78075,0.0401882,-8.76075e-06,-1.97193e-08,1.12783e-11,81164.5,17.38], Tmin=(100,'K'), Tmax=(955.832,'K')), NASAPolynomial(coeffs=[12.0562,0.0178508,-6.13458e-06,1.06669e-09,-7.39864e-14,78256.3,-36.6668], Tmin=(955.832,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(674.111,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(224.491,'J/(mol*K)'), comment="""Thermo group additivity estimation: group(Cs-(Cds-Cds)HHH) + group(Cds-Cds(Cds-Cds)Cs) + group(Cds-Cds(Cds-Cds)H) + group(Cds-CdsHH) + group(Cds-CdsHH) + radical(Cds_P) + radical(AllylJ2_triplet)"""),
)

species(
    label = '[CH][C]=[CH](21256)',
    structure = SMILES('[CH][C]=[CH]'),
    E0 = (861.746,'kJ/mol'),
    modes = [
        HarmonicOscillator(frequencies=([1685,370,3120,650,792.5,1650,180,180],'cm^-1')),
        HinderedRotor(inertia=(2.1891,'amu*angstrom^2'), symmetry=1, barrier=(50.3317,'kJ/mol'), semiclassical=False),
    ],
    spinMultiplicity = 5,
    opticalIsomers = 1,
    molecularWeight = (38.048,'amu'),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[3.18317,0.0164338,-7.13252e-06,1.19383e-09,-3.27944e-14,103675,12.0918], Tmin=(100,'K'), Tmax=(1799.19,'K')), NASAPolynomial(coeffs=[6.32962,0.0112581,-4.33439e-06,7.19107e-10,-4.49321e-14,102248,-5.75439], Tmin=(1799.19,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(861.746,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(103.931,'J/(mol*K)'), comment="""Thermo group additivity estimation: group(Cs-(Cds-Cds)HHH) + group(Cds-CdsCsH) + group(Cds-CdsHH) + radical(Cds_P) + radical(Cds_S) + radical(AllylJ2_triplet)"""),
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
    label = '[CH]C([CH])=C([CH2])C(19692)',
    structure = SMILES('[CH]C([CH])=C([CH2])C'),
    E0 = (784.856,'kJ/mol'),
    modes = [
        HarmonicOscillator(frequencies=([2750,2800,2850,1350,1500,750,1050,1375,1000,325,375,415,465,420,450,1700,1750,3000,3100,440,815,1455,1000,200,800,1000,1200,1400,1600],'cm^-1')),
        HinderedRotor(inertia=(0.156089,'amu*angstrom^2'), symmetry=1, barrier=(3.5888,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.156089,'amu*angstrom^2'), symmetry=1, barrier=(3.5888,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.156089,'amu*angstrom^2'), symmetry=1, barrier=(3.5888,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.156089,'amu*angstrom^2'), symmetry=1, barrier=(3.5888,'kJ/mol'), semiclassical=False),
    ],
    spinMultiplicity = 6,
    opticalIsomers = 1,
    molecularWeight = (79.1198,'amu'),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[0.948164,0.0613832,-3.83155e-05,1.23014e-08,-1.65031e-12,94510.7,23.945], Tmin=(100,'K'), Tmax=(1662.91,'K')), NASAPolynomial(coeffs=[12.3948,0.0338493,-1.3479e-05,2.34436e-09,-1.53386e-13,90703.8,-37.0999], Tmin=(1662.91,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(784.856,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(291.007,'J/(mol*K)'), comment="""Thermo group additivity estimation: group(Cs-(Cds-Cds)HHH) + group(Cs-(Cds-Cds)HHH) + group(Cs-(Cds-Cds)HHH) + group(Cs-(Cds-Cds)HHH) + group(Cds-CdsCsCs) + group(Cds-CdsCsCs) + radical(AllylJ2_triplet) + radical(AllylJ2_triplet) + radical(Allyl_P)"""),
)

species(
    label = '[CH]C([CH2])=C([CH2])[CH2](19199)',
    structure = SMILES('[CH]C([CH2])=C([CH2])[CH2]'),
    E0 = (717.17,'kJ/mol'),
    modes = [
        HarmonicOscillator(frequencies=([325,375,415,465,420,450,1700,1750,3000,3020,3040,3060,3080,3100,415,440,465,780,815,850,1435,1455,1475,900,1000,1100,457.373,457.42,457.49],'cm^-1')),
        HinderedRotor(inertia=(0.34168,'amu*angstrom^2'), symmetry=1, barrier=(50.7278,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.341455,'amu*angstrom^2'), symmetry=1, barrier=(50.7291,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.341283,'amu*angstrom^2'), symmetry=1, barrier=(50.7299,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.341662,'amu*angstrom^2'), symmetry=1, barrier=(50.7285,'kJ/mol'), semiclassical=False),
    ],
    spinMultiplicity = 6,
    opticalIsomers = 1,
    molecularWeight = (79.1198,'amu'),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[0.976316,0.056264,-2.15955e-05,-1.26555e-08,9.05162e-12,86373.5,23.2889], Tmin=(100,'K'), Tmax=(998.136,'K')), NASAPolynomial(coeffs=[13.9159,0.0264994,-1.00626e-05,1.81444e-09,-1.26205e-13,82690,-44.6253], Tmin=(998.136,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(717.17,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(291.007,'J/(mol*K)'), comment="""Thermo group additivity estimation: group(Cs-(Cds-Cds)HHH) + group(Cs-(Cds-Cds)HHH) + group(Cs-(Cds-Cds)HHH) + group(Cs-(Cds-Cds)HHH) + group(Cds-CdsCsCs) + group(Cds-CdsCsCs) + radical(Allyl_P) + radical(Allyl_P) + radical(AllylJ2_triplet) + radical(Allyl_P)"""),
)

species(
    label = '[CH]C([CH])=C([CH2])[CH2](21436)',
    structure = SMILES('[CH]C([CH])=C([CH2])[CH2]'),
    E0 = (936.355,'kJ/mol'),
    modes = [
        HarmonicOscillator(frequencies=([325,375,415,465,420,450,1700,1750,3000,3033.33,3066.67,3100,415,465,780,850,1435,1475,900,1100,200,800,1000,1200,1400,1600],'cm^-1')),
        HinderedRotor(inertia=(0.156089,'amu*angstrom^2'), symmetry=1, barrier=(3.5888,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.156089,'amu*angstrom^2'), symmetry=1, barrier=(3.5888,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.156089,'amu*angstrom^2'), symmetry=1, barrier=(3.5888,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.156089,'amu*angstrom^2'), symmetry=1, barrier=(3.5888,'kJ/mol'), semiclassical=False),
    ],
    spinMultiplicity = 7,
    opticalIsomers = 1,
    molecularWeight = (78.1118,'amu'),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[0.928701,0.0592111,-3.14168e-05,1.94434e-09,2.68939e-12,112735,24.29], Tmin=(100,'K'), Tmax=(1076.6,'K')), NASAPolynomial(coeffs=[11.9821,0.0319589,-1.26959e-05,2.27125e-09,-1.54377e-13,109554,-33.5706], Tmin=(1076.6,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(936.355,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(266.063,'J/(mol*K)'), comment="""Thermo group additivity estimation: group(Cs-(Cds-Cds)HHH) + group(Cs-(Cds-Cds)HHH) + group(Cs-(Cds-Cds)HHH) + group(Cs-(Cds-Cds)HHH) + group(Cds-CdsCsCs) + group(Cds-CdsCsCs) + radical(Allyl_P) + radical(AllylJ2_triplet) + radical(AllylJ2_triplet) + radical(Allyl_P)"""),
)

species(
    label = '[CH]C(=[CH])C(=C)C(19687)',
    structure = SMILES('[CH]C(=[CH])C(=C)C'),
    E0 = (636.521,'kJ/mol'),
    modes = [
        HarmonicOscillator(frequencies=([3120,650,792.5,1650,2750,2800,2850,1350,1500,750,1050,1375,1000,325,375,415,465,420,450,1700,1750,2950,3100,1380,975,1025,1650,180,180,180],'cm^-1')),
        HinderedRotor(inertia=(2.11706,'amu*angstrom^2'), symmetry=1, barrier=(48.6754,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(2.11657,'amu*angstrom^2'), symmetry=1, barrier=(48.6641,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(2.11677,'amu*angstrom^2'), symmetry=1, barrier=(48.6687,'kJ/mol'), semiclassical=False),
    ],
    spinMultiplicity = 4,
    opticalIsomers = 1,
    molecularWeight = (79.1198,'amu'),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[0.895327,0.0606763,-4.09674e-05,9.40654e-09,1.0466e-12,76674.1,21.159], Tmin=(100,'K'), Tmax=(1031.13,'K')), NASAPolynomial(coeffs=[13.2617,0.0259733,-9.78747e-06,1.72745e-09,-1.17364e-13,73418.4,-42.3016], Tmin=(1031.13,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(636.521,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(295.164,'J/(mol*K)'), comment="""Thermo group additivity estimation: group(Cs-(Cds-Cds)HHH) + group(Cs-(Cds-Cds)HHH) + group(Cds-Cds(Cds-Cds)Cs) + group(Cds-Cds(Cds-Cds)Cs) + group(Cds-CdsHH) + group(Cds-CdsHH) + radical(AllylJ2_triplet) + radical(Cds_P)"""),
)

species(
    label = '[CH]C(=C)C([CH2])=C(18076)',
    structure = SMILES('[CH]C(=C)C([CH2])=C'),
    E0 = (540.924,'kJ/mol'),
    modes = [
        HarmonicOscillator(frequencies=([2950,3000,3050,3100,1330,1430,900,1050,1000,1050,1600,1700,325,375,415,465,420,450,1700,1750,3000,3100,440,815,1455,1000,195.967,195.967,195.967,195.967],'cm^-1')),
        HinderedRotor(inertia=(1.80854,'amu*angstrom^2'), symmetry=1, barrier=(49.2857,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(1.80854,'amu*angstrom^2'), symmetry=1, barrier=(49.2857,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(1.80854,'amu*angstrom^2'), symmetry=1, barrier=(49.2857,'kJ/mol'), semiclassical=False),
    ],
    spinMultiplicity = 4,
    opticalIsomers = 1,
    molecularWeight = (79.1198,'amu'),
    collisionModel = TransportData(shapeIndex=2, epsilon=(3249.64,'J/mol'), sigma=(5.69809,'angstroms'), dipoleMoment=(0,'C*m'), polarizability=(0,'angstroms^3'), rotrelaxcollnum=0, comment="""Epsilon & sigma estimated with Tc=507.59 K, Pc=39.86 bar (from Joback method)"""),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[1.01928,0.0540671,-1.36977e-05,-2.39883e-08,1.40947e-11,65175.7,20.5913], Tmin=(100,'K'), Tmax=(963.045,'K')), NASAPolynomial(coeffs=[14.94,0.0233813,-8.16493e-06,1.43762e-09,-1.00384e-13,61236.2,-52.5766], Tmin=(963.045,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(540.924,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(295.164,'J/(mol*K)'), comment="""Thermo group additivity estimation: group(Cs-(Cds-Cds)HHH) + group(Cs-(Cds-Cds)HHH) + group(Cds-Cds(Cds-Cds)Cs) + group(Cds-Cds(Cds-Cds)Cs) + group(Cds-CdsHH) + group(Cds-CdsHH) + radical(Allyl_P) + radical(AllylJ2_triplet)"""),
)

species(
    label = '[CH]C(=[CH])C[CH][CH2](19224)',
    structure = SMILES('[CH]C(=[CH])C[CH][CH2]'),
    E0 = (938.116,'kJ/mol'),
    modes = [
        HarmonicOscillator(frequencies=([3000,3100,440,815,1455,1000,2750,2850,1437.5,1250,1305,750,350,350,440,435,1725,3120,650,792.5,1650,3025,407.5,1350,352.5,615.299,615.3,615.305,615.312],'cm^-1')),
        HinderedRotor(inertia=(0.000445249,'amu*angstrom^2'), symmetry=1, barrier=(0.119627,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.202719,'amu*angstrom^2'), symmetry=1, barrier=(54.4629,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.202717,'amu*angstrom^2'), symmetry=1, barrier=(54.4628,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.202714,'amu*angstrom^2'), symmetry=1, barrier=(54.4629,'kJ/mol'), semiclassical=False),
    ],
    spinMultiplicity = 6,
    opticalIsomers = 1,
    molecularWeight = (79.1198,'amu'),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[1.60609,0.0559424,-4.67467e-05,2.62878e-08,-6.96801e-12,112913,27.0253], Tmin=(100,'K'), Tmax=(855.646,'K')), NASAPolynomial(coeffs=[5.22312,0.039034,-1.71061e-05,3.19437e-09,-2.20863e-13,112294,10.1391], Tmin=(855.646,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(938.116,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(291.007,'J/(mol*K)'), comment="""Thermo group additivity estimation: group(Cs-CsCsHH) + group(Cs-(Cds-Cds)CsHH) + group(Cs-CsHHH) + group(Cs-(Cds-Cds)HHH) + group(Cds-CdsCsCs) + group(Cds-CdsHH) + radical(RCCJC) + radical(AllylJ2_triplet) + radical(RCCJ) + radical(Cds_P)"""),
)

species(
    label = '[CH]C(=[CH])C1CC1(21437)',
    structure = SMILES('[CH]C(=[CH])C1CC1'),
    E0 = (693.048,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
    molecularWeight = (79.1198,'amu'),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[1.71395,0.034658,3.558e-05,-7.08242e-08,2.96133e-11,83450.8,21.4116], Tmin=(100,'K'), Tmax=(960.347,'K')), NASAPolynomial(coeffs=[13.4873,0.024205,-8.36016e-06,1.51583e-09,-1.09579e-13,79410.3,-44.1753], Tmin=(960.347,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(693.048,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(299.321,'J/(mol*K)'), comment="""Thermo group additivity estimation: group(Cs-(Cds-Cds)CsCsH) + group(Cs-CsCsHH) + group(Cs-CsCsHH) + group(Cs-(Cds-Cds)HHH) + group(Cds-CdsCsCs) + group(Cds-CdsHH) + ring(Cyclopropane) + radical(AllylJ2_triplet) + radical(Cds_P)"""),
)

species(
    label = '[CH]C1=CCC1[CH2](21438)',
    structure = SMILES('[CH]C1=CCC1[CH2]'),
    E0 = (649.427,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
    molecularWeight = (79.1198,'amu'),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[1.91665,0.0303916,4.52062e-05,-7.9824e-08,3.29846e-11,78197.1,21.7315], Tmin=(100,'K'), Tmax=(943.087,'K')), NASAPolynomial(coeffs=[12.3117,0.0250739,-8.00312e-06,1.3821e-09,-9.7868e-14,74512.2,-36.9507], Tmin=(943.087,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(649.427,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(299.321,'J/(mol*K)'), comment="""Thermo group additivity estimation: group(Cs-(Cds-Cds)CsCsH) + group(Cs-(Cds-Cds)CsHH) + group(Cs-CsHHH) + group(Cs-(Cds-Cds)HHH) + group(Cds-CdsCsCs) + group(Cds-CdsCsH) + ring(Cyclobutene) + radical(Isobutyl) + radical(AllylJ2_triplet)"""),
)

species(
    label = '[CH]C([CH])=C[CH2](21258)',
    structure = SMILES('[CH]C([CH])=C[CH2]'),
    E0 = (823.911,'kJ/mol'),
    modes = [
        HarmonicOscillator(frequencies=([3000,3100,440,815,1455,1000,3010,987.5,1337.5,450,1655,350,440,435,1725,328.03,328.033,328.034,328.035,328.036,328.04],'cm^-1')),
        HinderedRotor(inertia=(0.664758,'amu*angstrom^2'), symmetry=1, barrier=(50.762,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.664777,'amu*angstrom^2'), symmetry=1, barrier=(50.7619,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.664772,'amu*angstrom^2'), symmetry=1, barrier=(50.7618,'kJ/mol'), semiclassical=False),
    ],
    spinMultiplicity = 6,
    opticalIsomers = 1,
    molecularWeight = (65.0932,'amu'),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[1.78722,0.0415433,-8.93762e-06,-1.22359e-08,6.20359e-12,99179.3,20.1709], Tmin=(100,'K'), Tmax=(1059.47,'K')), NASAPolynomial(coeffs=[8.69747,0.0295269,-1.18496e-05,2.13401e-09,-1.45711e-13,96925.3,-17.2938], Tmin=(1059.47,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(823.911,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(220.334,'J/(mol*K)'), comment="""Thermo group additivity estimation: group(Cs-(Cds-Cds)HHH) + group(Cs-(Cds-Cds)HHH) + group(Cs-(Cds-Cds)HHH) + group(Cds-CdsCsCs) + group(Cds-CdsCsH) + radical(AllylJ2_triplet) + radical(Allyl_P) + radical(AllylJ2_triplet)"""),
)

species(
    label = '[CH]C(=[CH])C([CH])[CH2](19915)',
    structure = SMILES('[CH]C(=[CH])C([CH])[CH2]'),
    E0 = (1187.93,'kJ/mol'),
    modes = [
        HarmonicOscillator(frequencies=([3000,3100,440,815,1455,1000,3120,650,792.5,1650,1380,1390,370,380,2900,435,350,440,435,1725,200,800,1000,1200,1400,1600],'cm^-1')),
        HinderedRotor(inertia=(0.156089,'amu*angstrom^2'), symmetry=1, barrier=(3.5888,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.156089,'amu*angstrom^2'), symmetry=1, barrier=(3.5888,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.156089,'amu*angstrom^2'), symmetry=1, barrier=(3.5888,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.156089,'amu*angstrom^2'), symmetry=1, barrier=(3.5888,'kJ/mol'), semiclassical=False),
    ],
    spinMultiplicity = 7,
    opticalIsomers = 1,
    molecularWeight = (78.1118,'amu'),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[0.940463,0.0619128,-5.53295e-05,2.71444e-08,-5.40431e-12,142990,27.4197], Tmin=(100,'K'), Tmax=(1209.48,'K')), NASAPolynomial(coeffs=[12.3041,0.0243309,-8.7201e-06,1.45308e-09,-9.3883e-14,140241,-29.5642], Tmin=(1209.48,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(1187.93,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(266.063,'J/(mol*K)'), comment="""Thermo group additivity estimation: group(Cs-(Cds-Cds)CsCsH) + group(Cs-CsHHH) + group(Cs-CsHHH) + group(Cs-(Cds-Cds)HHH) + group(Cds-CdsCsCs) + longDistanceInteraction_noncyclic(CdCs-ST) + group(Cds-CdsHH) + radical(CCJ2_triplet) + radical(Isobutyl) + radical(Cds_P) + radical(AllylJ2_triplet)"""),
)

species(
    label = '[C]C(=[CH])C([CH2])[CH2](21439)',
    structure = SMILES('[C]C(=[CH])C([CH2])[CH2]'),
    E0 = (1243.59,'kJ/mol'),
    modes = [
        HarmonicOscillator(frequencies=([3000,3033.33,3066.67,3100,415,465,780,850,1435,1475,900,1100,3120,650,792.5,1650,1380,1390,370,380,2900,435,350,440,435,1725,367.708],'cm^-1')),
        HinderedRotor(inertia=(0.00122502,'amu*angstrom^2'), symmetry=1, barrier=(0.119627,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.0701364,'amu*angstrom^2'), symmetry=1, barrier=(6.83488,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.706325,'amu*angstrom^2'), symmetry=1, barrier=(68.0055,'kJ/mol'), semiclassical=False),
    ],
    spinMultiplicity = 7,
    opticalIsomers = 1,
    molecularWeight = (78.1118,'amu'),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[1.15357,0.0603236,-5.98416e-05,2.62891e-08,-2.27683e-12,149675,24.7839], Tmin=(100,'K'), Tmax=(839.269,'K')), NASAPolynomial(coeffs=[12.9498,0.0168364,-4.87811e-06,7.08259e-10,-4.21914e-14,147246,-32.7311], Tmin=(839.269,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(1243.59,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(270.22,'J/(mol*K)'), comment="""Thermo group additivity estimation: group(Cs-(Cds-Cds)CsCsH) + group(Cs-CsHHH) + group(Cs-CsHHH) + group(Cs-(Cds-Cds)HHH) + group(Cds-CdsCsCs) + longDistanceInteraction_noncyclic(CdCs-ST) + group(Cds-CdsHH) + radical(Cds_P) + radical(Isobutyl) + radical(Isobutyl) + radical(CJ3)"""),
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
    E0 = (944.8,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS2',
    E0 = (1011.11,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS3',
    E0 = (1128.03,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS4',
    E0 = (1065.41,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS5',
    E0 = (1039.88,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS6',
    E0 = (998.658,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS7',
    E0 = (1086.29,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS8',
    E0 = (1136.93,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS9',
    E0 = (1346.04,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS10',
    E0 = (1148.16,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS11',
    E0 = (967.662,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS12',
    E0 = (1008.2,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS13',
    E0 = (1102.12,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS14',
    E0 = (952.889,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS15',
    E0 = (953.085,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS16',
    E0 = (1239.6,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS17',
    E0 = (1399.74,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS18',
    E0 = (1455.4,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

reaction(
    label = 'reaction1',
    reactants = ['[CH]C(=[CH])C([CH2])[CH2](19200)'],
    products = ['[CH2]C=C(87)', '[CH]=C=[CH](18734)'],
    transitionState = 'TS1',
    kinetics = Arrhenius(A=(5e+12,'s^-1'), n=0, Ea=(0,'kJ/mol'), T0=(1,'K'), Tmin=(300,'K'), Tmax=(1500,'K'), comment="""Exact match found for rate rule [RJJ]
Euclidian distance = 0
family: 1,4_Linear_birad_scission"""),
)

reaction(
    label = 'reaction2',
    reactants = ['H(8)', '[CH]C(=[CH])C([CH2])=C(19688)'],
    products = ['[CH]C(=[CH])C([CH2])[CH2](19200)'],
    transitionState = 'TS2',
    kinetics = Arrhenius(A=(170.395,'m^3/(mol*s)'), n=1.5621, Ea=(11.2886,'kJ/mol'), T0=(1,'K'), comment="""Estimated using an average for rate rule [Cds_Cds;HJ]
Euclidian distance = 0
family: R_Addition_MultipleBond"""),
)

reaction(
    label = 'reaction3',
    reactants = ['[CH](2815)', 'C#CC([CH2])[CH2](6419)'],
    products = ['[CH]C(=[CH])C([CH2])[CH2](19200)'],
    transitionState = 'TS3',
    kinetics = Arrhenius(A=(18.899,'m^3/(mol*s)'), n=1.76329, Ea=(16.1554,'kJ/mol'), T0=(1,'K'), comment="""Estimated using template [Ct-Cs_Ct-H;YJ] for rate rule [Ct-Cs_Ct-H;CH_quartet]
Euclidian distance = 2.0
family: R_Addition_MultipleBond"""),
)

reaction(
    label = 'reaction4',
    reactants = ['CH2(T)(28)', '[CH]C(=[CH])C=C(19261)'],
    products = ['[CH]C(=[CH])C([CH2])[CH2](19200)'],
    transitionState = 'TS4',
    kinetics = Arrhenius(A=(2.39248,'m^3/(mol*s)'), n=1.9454, Ea=(9.93359,'kJ/mol'), T0=(1,'K'), comment="""Estimated using template [Cds-OneDeH_Cds;YJ] for rate rule [Cds-OneDeH_Cds;CH2_triplet]
Euclidian distance = 2.0
family: R_Addition_MultipleBond"""),
)

reaction(
    label = 'reaction5',
    reactants = ['[CH2]C=C(87)', '[CH][C]=[CH](21256)'],
    products = ['[CH]C(=[CH])C([CH2])[CH2](19200)'],
    transitionState = 'TS5',
    kinetics = Arrhenius(A=(0.00179511,'m^3/(mol*s)'), n=2.50446, Ea=(20.5128,'kJ/mol'), T0=(1,'K'), comment="""Estimated using an average for rate rule [Cds_Cds;CJ]
Euclidian distance = 0
family: R_Addition_MultipleBond"""),
)

reaction(
    label = 'reaction6',
    reactants = ['[CH2][CH][CH2](6136)', '[CH]=C=[CH](18734)'],
    products = ['[CH]C(=[CH])C([CH2])[CH2](19200)'],
    transitionState = 'TS6',
    kinetics = Arrhenius(A=(0.523563,'m^3/(mol*s)'), n=2.10494, Ea=(22.6844,'kJ/mol'), T0=(1,'K'), comment="""Estimated using an average for rate rule [Ct_Ct;CJ]
Euclidian distance = 0
family: R_Addition_MultipleBond"""),
)

reaction(
    label = 'reaction7',
    reactants = ['[CH]C(=[CH])C([CH2])[CH2](19200)'],
    products = ['[CH]C([CH])=C([CH2])C(19692)'],
    transitionState = 'TS7',
    kinetics = Arrhenius(A=(34962.3,'s^-1'), n=2.56136, Ea=(141.488,'kJ/mol'), T0=(1,'K'), comment="""Estimated using an average for rate rule [R2H_S;C_rad_out_2H;XH_out]
Euclidian distance = 0
Multiplied by reaction path degeneracy 2.0
family: intra_H_migration"""),
)

reaction(
    label = 'reaction8',
    reactants = ['[CH]C(=[CH])C([CH2])[CH2](19200)'],
    products = ['[CH]C([CH2])=C([CH2])[CH2](19199)'],
    transitionState = 'TS8',
    kinetics = Arrhenius(A=(13437.7,'s^-1'), n=2.58467, Ea=(192.129,'kJ/mol'), T0=(1,'K'), comment="""Estimated using an average for rate rule [R3H_DS;Cd_rad_out_singleH;XH_out]
Euclidian distance = 0
Multiplied by reaction path degeneracy 2.0
family: intra_H_migration"""),
)

reaction(
    label = 'reaction9',
    reactants = ['[CH2][CH][CH2](6136)', '[CH][C]=[CH](21256)'],
    products = ['[CH]C(=[CH])C([CH2])[CH2](19200)'],
    transitionState = 'TS9',
    kinetics = Arrhenius(A=(1.9789e+07,'m^3/(mol*s)'), n=-0.126319, Ea=(0,'kJ/mol'), T0=(1,'K'), comment="""Estimated using an average for rate rule [Y_rad;Y_rad]
Euclidian distance = 0
family: R_Recombination
Ea raised from -15.6 to -15.6 kJ/mol.
Ea raised from -15.6 to 0 kJ/mol."""),
)

reaction(
    label = 'reaction10',
    reactants = ['H(8)', '[CH]C([CH])=C([CH2])[CH2](21436)'],
    products = ['[CH]C(=[CH])C([CH2])[CH2](19200)'],
    transitionState = 'TS10',
    kinetics = Arrhenius(A=(4.34078e+06,'m^3/(mol*s)'), n=0.278577, Ea=(0,'kJ/mol'), T0=(1,'K'), comment="""Estimated using an average for rate rule [Y_rad;H_rad]
Euclidian distance = 0
family: R_Recombination
Ea raised from -1.4 to 0 kJ/mol."""),
)

reaction(
    label = 'reaction11',
    reactants = ['[CH]C(=[CH])C([CH2])[CH2](19200)'],
    products = ['[CH]C(=[CH])C(=C)C(19687)'],
    transitionState = 'TS11',
    kinetics = Arrhenius(A=(3.898e+11,'s^-1'), n=0.486, Ea=(22.8614,'kJ/mol'), T0=(1,'K'), comment="""Estimated using an average for rate rule [R2radExo;Y_rad;XH_Rrad]
Euclidian distance = 0
Multiplied by reaction path degeneracy 2.0
family: Intra_Disproportionation"""),
)

reaction(
    label = 'reaction12',
    reactants = ['[CH]C(=[CH])C([CH2])[CH2](19200)'],
    products = ['[CH]C(=C)C([CH2])=C(18076)'],
    transitionState = 'TS12',
    kinetics = Arrhenius(A=(1.4874e+09,'s^-1'), n=1.045, Ea=(63.4002,'kJ/mol'), T0=(1,'K'), comment="""Estimated using an average for rate rule [R3radExo;Y_rad;XH_Rrad]
Euclidian distance = 0
Multiplied by reaction path degeneracy 2.0
family: Intra_Disproportionation"""),
)

reaction(
    label = 'reaction13',
    reactants = ['[CH]C(=[CH])C([CH2])[CH2](19200)'],
    products = ['[CH]C(=[CH])C[CH][CH2](19224)'],
    transitionState = 'TS13',
    kinetics = Arrhenius(A=(2.66e+08,'s^-1'), n=1.36, Ea=(157.318,'kJ/mol'), T0=(1,'K'), comment="""Estimated using an average for rate rule [cCsCJ;CsJ-HH;C]
Euclidian distance = 0
Multiplied by reaction path degeneracy 2.0
family: 1,2_shiftC"""),
)

reaction(
    label = 'reaction14',
    reactants = ['[CH]C(=[CH])C([CH2])[CH2](19200)'],
    products = ['[CH]C(=[CH])C1CC1(21437)'],
    transitionState = 'TS14',
    kinetics = Arrhenius(A=(7.38971e+10,'s^-1'), n=0.0476667, Ea=(8.08907,'kJ/mol'), T0=(1,'K'), comment="""Estimated using template [Rn;C_rad_out_2H;Cpri_rad_out_2H] for rate rule [R3_SS;C_rad_out_2H;Cpri_rad_out_2H]
Euclidian distance = 2.0
family: Birad_recombination"""),
)

reaction(
    label = 'reaction15',
    reactants = ['[CH]C(=[CH])C([CH2])[CH2](19200)'],
    products = ['[CH]C1=CCC1[CH2](21438)'],
    transitionState = 'TS15',
    kinetics = Arrhenius(A=(6.48e+12,'s^-1'), n=-0.305, Ea=(8.28432,'kJ/mol'), T0=(1,'K'), comment="""Estimated using template [R4;C_rad_out_2H;Ypri_rad_out] for rate rule [R4_SSD;C_rad_out_2H;CdsinglepriH_rad_out]
Euclidian distance = 2.2360679775
Multiplied by reaction path degeneracy 4.0
family: Birad_recombination"""),
)

reaction(
    label = 'reaction16',
    reactants = ['CH2(T)(28)', '[CH]C([CH])=C[CH2](21258)'],
    products = ['[CH]C(=[CH])C([CH2])[CH2](19200)'],
    transitionState = 'TS16',
    kinetics = Arrhenius(A=(1.14854e+06,'m^3/(mol*s)'), n=0.575199, Ea=(34.3157,'kJ/mol'), T0=(1,'K'), comment="""Estimated using an average for rate rule [Y_rad;Birad]
Euclidian distance = 0
family: Birad_R_Recombination"""),
)

reaction(
    label = 'reaction17',
    reactants = ['H(8)', '[CH]C(=[CH])C([CH])[CH2](19915)'],
    products = ['[CH]C(=[CH])C([CH2])[CH2](19200)'],
    transitionState = 'TS17',
    kinetics = Arrhenius(A=(1e+07,'m^3/(mol*s)'), n=0, Ea=(0,'kJ/mol'), T0=(1,'K'), comment="""Estimated using an average for rate rule [H_rad;Birad]
Euclidian distance = 0
family: Birad_R_Recombination"""),
)

reaction(
    label = 'reaction18',
    reactants = ['H(8)', '[C]C(=[CH])C([CH2])[CH2](21439)'],
    products = ['[CH]C(=[CH])C([CH2])[CH2](19200)'],
    transitionState = 'TS18',
    kinetics = Arrhenius(A=(1e+07,'m^3/(mol*s)'), n=0, Ea=(0,'kJ/mol'), T0=(1,'K'), comment="""Estimated using an average for rate rule [H_rad;Birad]
Euclidian distance = 0
family: Birad_R_Recombination"""),
)

network(
    label = '4465',
    isomers = [
        '[CH]C(=[CH])C([CH2])[CH2](19200)',
    ],
    reactants = [
        ('[CH2]C=C(87)', '[CH]=C=[CH](18734)'),
    ],
    bathGas = {
        'N2': 0.25,
        'Ne': 0.25,
        'He': 0.25,
        'Ar': 0.25,
    },
)

pressureDependence(
    label = '4465',
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

