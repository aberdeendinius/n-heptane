species(
    label = '[CH]=[C]C=C([CH2])[CH]O(21903)',
    structure = SMILES('[CH]=[C]C=C([CH2])[CH]O'),
    E0 = (547.147,'kJ/mol'),
    modes = [
        HarmonicOscillator(frequencies=([3120,650,792.5,1650,3615,1277.5,1000,3025,407.5,1350,352.5,350,440,435,1725,3000,3100,440,815,1455,1000,1685,370,3010,987.5,1337.5,450,1655,433.307],'cm^-1')),
        HinderedRotor(inertia=(0.50806,'amu*angstrom^2'), symmetry=1, barrier=(67.0981,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.507427,'amu*angstrom^2'), symmetry=1, barrier=(67.0725,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.224943,'amu*angstrom^2'), symmetry=1, barrier=(29.9272,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.501106,'amu*angstrom^2'), symmetry=1, barrier=(67.0762,'kJ/mol'), semiclassical=False),
    ],
    spinMultiplicity = 5,
    opticalIsomers = 1,
    molecularWeight = (94.1112,'amu'),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[0.902672,0.0580589,-3.00243e-05,-1.15182e-08,1.15877e-11,65927.4,27.5855], Tmin=(100,'K'), Tmax=(922.595,'K')), NASAPolynomial(coeffs=[16.6204,0.0160842,-4.32983e-06,6.61554e-10,-4.43179e-14,61913.4,-53.0132], Tmin=(922.595,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(547.147,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(291.007,'J/(mol*K)'), comment="""Thermo group additivity estimation: group(O2s-CsH) + group(Cs-(Cds-Cds)OsHH) + group(Cs-(Cds-Cds)HHH) + group(Cds-CdsCsCs) + group(Cds-Cds(Cds-Cds)H) + group(Cds-Cds(Cds-Cds)H) + group(Cds-CdsHH) + radical(Cds_P) + radical(C=CCJO) + radical(C=CC=CCJ) + radical(C=CJC=C)"""),
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
    label = '[CH2][C]([CH]O)C1[C]=C1(29830)',
    structure = SMILES('[CH2][C]([CH]O)C1[C]=C1'),
    E0 = (790.825,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
    molecularWeight = (94.1112,'amu'),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[0.67854,0.0764542,-0.000105401,8.03814e-08,-2.43354e-11,95230.6,26.483], Tmin=(100,'K'), Tmax=(878.224,'K')), NASAPolynomial(coeffs=[10.6559,0.0259704,-1.05665e-05,1.85674e-09,-1.21898e-13,93672.5,-19.25], Tmin=(878.224,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(790.825,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(291.007,'J/(mol*K)'), comment="""Thermo group additivity estimation: group(O2s-CsH) + group(Cs-CsCsCsH) + group(Cs-(Cds-Cds)(Cds-Cds)CsH) + group(Cs-CsOsHH) + group(Cs-CsHHH) + group(Cds-CdsCsH) + group(Cds-CdsCsH) + ring(Cyclopropene) + radical(CCsJOH) + radical(cyclopropenyl-vinyl) + radical(CCJ(C)CO) + radical(Isobutyl)"""),
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
    label = '[CH]=C=CC([CH2])=C[O](21884)',
    structure = SMILES('[CH]=C=CC([CH2])=C[O]'),
    E0 = (435.555,'kJ/mol'),
    modes = [
        HarmonicOscillator(frequencies=([3000,3100,440,815,1455,1000,3120,650,792.5,1650,540,610,2055,350,440,435,1725,2995,3025,975,1000,1300,1375,400,500,1630,1680,180],'cm^-1')),
        HinderedRotor(inertia=(1.55659,'amu*angstrom^2'), symmetry=1, barrier=(35.7891,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(1.55719,'amu*angstrom^2'), symmetry=1, barrier=(35.803,'kJ/mol'), semiclassical=False),
    ],
    spinMultiplicity = 4,
    opticalIsomers = 1,
    molecularWeight = (93.1033,'amu'),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[0.626417,0.0593098,-2.25474e-05,-3.25963e-08,2.28341e-11,52520.7,23.7429], Tmin=(100,'K'), Tmax=(902.146,'K')), NASAPolynomial(coeffs=[22.4044,0.00362593,2.07193e-06,-5.63809e-10,3.90027e-14,46927.8,-88.3], Tmin=(902.146,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(435.555,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(274.378,'J/(mol*K)'), comment="""Thermo group additivity estimation: group(O2s-(Cds-Cd)H) + group(Cs-(Cds-Cds)HHH) + group(Cds-Cds(Cds-Cds)Cs) + group(Cds-Cds(Cds-Cds)H) + group(Cds-CdsOsH) + group(Cds-CdsHH) + group(Cdd-CdsCds) + radical(C=C=CJ) + radical(C=COJ) + radical(Allyl_P)"""),
)

species(
    label = '[CH]=C=C=C([CH2])[CH]O(29831)',
    structure = SMILES('[CH]=C=C=C([CH2])[CH]O'),
    E0 = (518.361,'kJ/mol'),
    modes = [
        HarmonicOscillator(frequencies=([3120,650,792.5,1650,3615,1277.5,1000,3025,407.5,1350,352.5,540,563.333,586.667,610,1970,2140,350,440,435,1725,3000,3100,440,815,1455,1000],'cm^-1')),
        HinderedRotor(inertia=(1.65823,'amu*angstrom^2'), symmetry=1, barrier=(38.1259,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(1.66127,'amu*angstrom^2'), symmetry=1, barrier=(38.1958,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(1.65616,'amu*angstrom^2'), symmetry=1, barrier=(38.0784,'kJ/mol'), semiclassical=False),
    ],
    spinMultiplicity = 4,
    opticalIsomers = 1,
    molecularWeight = (93.1033,'amu'),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[0.401001,0.0667925,-6.69682e-05,3.32717e-08,-6.31755e-12,62484.6,27.8232], Tmin=(100,'K'), Tmax=(1411.32,'K')), NASAPolynomial(coeffs=[17.8761,0.0119519,-3.03559e-06,4.04748e-10,-2.3097e-14,58081.1,-60.6302], Tmin=(1411.32,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(518.361,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(270.22,'J/(mol*K)'), comment="""Thermo group additivity estimation: group(O2s-CsH) + group(Cs-(Cds-Cds)OsHH) + group(Cs-(Cds-Cds)HHH) + group(Cds-CdsCsCs) + group(Cds-CdsHH) + group(Cdd-CdsCds) + group(Cdd-CdsCds) + radical(C=CCJO) + radical(C=C=CJ) + radical(Allyl_P)"""),
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
    label = '[CH]=C=C[C]=CO(23383)',
    structure = SMILES('[CH]=C=C[C]=CO'),
    E0 = (379.179,'kJ/mol'),
    modes = [
        HarmonicOscillator(frequencies=([3120,650,792.5,1650,1685,370,540,610,2055,3615,1277.5,1000,2995,3025,975,1000,1300,1375,400,500,1630,1680],'cm^-1')),
        HinderedRotor(inertia=(1.52208,'amu*angstrom^2'), symmetry=1, barrier=(34.9955,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(1.52206,'amu*angstrom^2'), symmetry=1, barrier=(34.9951,'kJ/mol'), semiclassical=False),
    ],
    spinMultiplicity = 3,
    opticalIsomers = 1,
    molecularWeight = (80.0847,'amu'),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[-0.147627,0.0687227,-8.0576e-05,4.32471e-08,-8.29888e-12,45773.5,23.8199], Tmin=(100,'K'), Tmax=(1567.29,'K')), NASAPolynomial(coeffs=[19.008,-0.000223152,4.60602e-06,-1.1515e-09,8.56534e-14,42232.5,-69.3436], Tmin=(1567.29,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(379.179,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(224.491,'J/(mol*K)'), comment="""Thermo group additivity estimation: group(O2s-(Cds-Cd)H) + group(Cds-Cds(Cds-Cds)H) + group(Cds-Cds(Cds-Cds)H) + group(Cds-CdsOsH) + group(Cds-CdsHH) + group(Cdd-CdsCds) + radical(C=CJC=C) + radical(C=C=CJ)"""),
)

species(
    label = '[CH]O(5471)',
    structure = SMILES('[CH]O'),
    E0 = (205.907,'kJ/mol'),
    modes = [
        HarmonicOscillator(frequencies=([3615,1277.5,1000,402.686,3356.18],'cm^-1')),
        HinderedRotor(inertia=(0.0105042,'amu*angstrom^2'), symmetry=1, barrier=(23.1306,'kJ/mol'), semiclassical=False),
    ],
    spinMultiplicity = 3,
    opticalIsomers = 1,
    molecularWeight = (30.026,'amu'),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[3.76003,0.0029575,8.86344e-06,-1.3392e-08,5.33433e-12,24775.7,6.76105], Tmin=(100,'K'), Tmax=(943.117,'K')), NASAPolynomial(coeffs=[5.07489,0.00326005,-9.68482e-07,1.67779e-10,-1.21779e-14,24266.2,-0.891576], Tmin=(943.117,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(205.907,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(78.9875,'J/(mol*K)'), comment="""Thermo group additivity estimation: group(O2s-CsH) + group(Cs-OsHHH) + radical(OsCsJ2H_triplet)"""),
)

species(
    label = '[CH]=C=C[C]=C(19277)',
    structure = SMILES('[CH]=C=C[C]=C'),
    E0 = (587.972,'kJ/mol'),
    modes = [
        HarmonicOscillator(frequencies=([3120,650,792.5,1650,2950,3100,1380,975,1025,1650,540,610,2055,3010,987.5,1337.5,450,1655,1685,370],'cm^-1')),
        HinderedRotor(inertia=(1.7607,'amu*angstrom^2'), symmetry=1, barrier=(40.4819,'kJ/mol'), semiclassical=False),
    ],
    spinMultiplicity = 3,
    opticalIsomers = 1,
    molecularWeight = (64.0853,'amu'),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[1.63528,0.044168,-4.52114e-05,2.40401e-08,-4.84455e-12,70808.6,17.5611], Tmin=(100,'K'), Tmax=(1402.09,'K')), NASAPolynomial(coeffs=[11.4691,0.0097722,-1.62986e-06,9.23713e-11,5.78346e-16,68674.3,-30.9823], Tmin=(1402.09,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(587.972,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(203.705,'J/(mol*K)'), comment="""Thermo group additivity estimation: group(Cds-Cds(Cds-Cds)H) + group(Cds-Cds(Cds-Cds)H) + group(Cds-CdsHH) + group(Cds-CdsHH) + group(Cdd-CdsCds) + radical(C=CJC=C) + radical(C=C=CJ)"""),
)

species(
    label = '[CH]=[C]C=C([CH2])C[O](22137)',
    structure = SMILES('[CH]=[C]C=C([CH2])C[O]'),
    E0 = (655.557,'kJ/mol'),
    modes = [
        HarmonicOscillator(frequencies=([3120,650,792.5,1650,350,440,435,1725,2750,2850,1437.5,1250,1305,750,350,3000,3100,440,815,1455,1000,1685,370,3010,987.5,1337.5,450,1655,180,470.194],'cm^-1')),
        HinderedRotor(inertia=(0.0149028,'amu*angstrom^2'), symmetry=1, barrier=(69.1341,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(3.00692,'amu*angstrom^2'), symmetry=1, barrier=(69.135,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(3.00856,'amu*angstrom^2'), symmetry=1, barrier=(69.1726,'kJ/mol'), semiclassical=False),
    ],
    spinMultiplicity = 5,
    opticalIsomers = 1,
    molecularWeight = (94.1112,'amu'),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[1.2047,0.0645813,-7.58703e-05,5.56196e-08,-1.71385e-11,78943.1,26.7682], Tmin=(100,'K'), Tmax=(817.437,'K')), NASAPolynomial(coeffs=[7.51007,0.0321696,-1.35367e-05,2.45233e-09,-1.65384e-13,77964.2,-2.06223], Tmin=(817.437,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(655.557,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(295.164,'J/(mol*K)'), comment="""Thermo group additivity estimation: group(O2s-CsH) + group(Cs-(Cds-Cds)OsHH) + group(Cs-(Cds-Cds)HHH) + group(Cds-CdsCsCs) + group(Cds-Cds(Cds-Cds)H) + group(Cds-Cds(Cds-Cds)H) + group(Cds-CdsHH) + radical(C=CJC=C) + radical(Cds_P) + radical(CCOJ) + radical(C=CC=CCJ)"""),
)

species(
    label = '[CH]C=[C]C([CH2])=CO(29832)',
    structure = SMILES('[CH]C=[C]C([CH2])=CO'),
    E0 = (532.691,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
    molecularWeight = (94.1112,'amu'),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[-1.13559,0.0862941,-8.65646e-05,4.19577e-08,-7.53389e-12,64275.8,29.6962], Tmin=(100,'K'), Tmax=(1613.39,'K')), NASAPolynomial(coeffs=[22.3663,0.0108611,-4.73659e-07,-2.10379e-10,2.20238e-14,58926.5,-88.0046], Tmin=(1613.39,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(532.691,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(291.007,'J/(mol*K)'), comment="""Thermo group additivity estimation: group(O2s-(Cds-Cd)H) + group(Cs-(Cds-Cds)HHH) + group(Cs-(Cds-Cds)HHH) + group(Cds-Cds(Cds-Cds)Cs) + group(Cds-CdsCsH) + group(Cds-Cds(Cds-Cds)H) + group(Cds-CdsOsH) + radical(AllylJ2_triplet) + radical(C=CJC=C) + radical(Allyl_P)"""),
)

species(
    label = '[CH]=[C][C]=C([CH2])CO(29833)',
    structure = SMILES('[CH]=[C][C]=C([CH2])CO'),
    E0 = (628.848,'kJ/mol'),
    modes = [
        HarmonicOscillator(frequencies=([3120,650,792.5,1650,3615,1277.5,1000,350,440,435,1725,2750,2850,1437.5,1250,1305,750,350,3000,3100,440,815,1455,1000,1670,1700,300,440,457.31],'cm^-1')),
        HinderedRotor(inertia=(0.392254,'amu*angstrom^2'), symmetry=1, barrier=(9.0187,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(2.77078,'amu*angstrom^2'), symmetry=1, barrier=(63.7057,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(2.76849,'amu*angstrom^2'), symmetry=1, barrier=(63.653,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(2.77103,'amu*angstrom^2'), symmetry=1, barrier=(63.7114,'kJ/mol'), semiclassical=False),
    ],
    spinMultiplicity = 5,
    opticalIsomers = 1,
    molecularWeight = (94.1112,'amu'),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[0.651168,0.0752798,-9.79207e-05,7.06528e-08,-2.00644e-11,75751.9,27.4861], Tmin=(100,'K'), Tmax=(955.065,'K')), NASAPolynomial(coeffs=[11.0884,0.0250929,-8.93106e-06,1.43793e-09,-8.87993e-14,74053.5,-20.842], Tmin=(955.065,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(628.848,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(291.007,'J/(mol*K)'), comment="""Thermo group additivity estimation: group(O2s-CsH) + group(Cs-(Cds-Cds)OsHH) + group(Cs-(Cds-Cds)HHH) + group(Cds-CdsCsCs) + group(Cds-Cds(Cds-Cds)H) + group(Cds-Cds(Cds-Cds)H) + group(Cds-CdsHH) + radical(Cds_P) + radical(C=CC=CCJ) + radical(C=CJC=C) + radical(C=CJC=C)"""),
)

species(
    label = '[CH][C]=[C]C(C)=CO(29834)',
    structure = SMILES('[CH][C]=[C]C(C)=CO'),
    E0 = (619.034,'kJ/mol'),
    modes = [
        HarmonicOscillator(frequencies=([1670,1700,300,440,2750,2800,2850,1350,1500,750,1050,1375,1000,3010,987.5,1337.5,450,1655,3615,1277.5,1000,350,440,435,1725,180,180,180,180],'cm^-1')),
        HinderedRotor(inertia=(2.01286,'amu*angstrom^2'), symmetry=1, barrier=(46.2797,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(2.01242,'amu*angstrom^2'), symmetry=1, barrier=(46.2694,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(2.01225,'amu*angstrom^2'), symmetry=1, barrier=(46.2656,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(2.01278,'amu*angstrom^2'), symmetry=1, barrier=(46.2778,'kJ/mol'), semiclassical=False),
    ],
    spinMultiplicity = 5,
    opticalIsomers = 1,
    molecularWeight = (94.1112,'amu'),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[-0.425212,0.0849329,-9.16428e-05,4.95841e-08,-1.02214e-11,74622.6,27.3187], Tmin=(100,'K'), Tmax=(1318.42,'K')), NASAPolynomial(coeffs=[19.9005,0.0152911,-3.33637e-06,3.43578e-10,-1.44165e-14,69956.1,-73.7313], Tmin=(1318.42,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(619.034,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(291.007,'J/(mol*K)'), comment="""Thermo group additivity estimation: group(O2s-(Cds-Cd)H) + group(Cs-(Cds-Cds)HHH) + group(Cs-(Cds-Cds)HHH) + group(Cds-Cds(Cds-Cds)Cs) + group(Cds-CdsCsH) + group(Cds-Cds(Cds-Cds)H) + group(Cds-CdsOsH) + radical(AllylJ2_triplet) + radical(C=CJC=C) + radical(Cds_S)"""),
)

species(
    label = '[CH2]C(=[C][C]=C)[CH]O(28964)',
    structure = SMILES('[CH2]C(=[C][C]=C)[CH]O'),
    E0 = (499.047,'kJ/mol'),
    modes = [
        HarmonicOscillator(frequencies=([3615,1277.5,1000,3025,407.5,1350,352.5,350,440,435,1725,2950,3100,1380,975,1025,1650,3000,3100,440,815,1455,1000,1670,1700,300,440,417.894,417.905],'cm^-1')),
        HinderedRotor(inertia=(0.447279,'amu*angstrom^2'), symmetry=1, barrier=(55.4633,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.447532,'amu*angstrom^2'), symmetry=1, barrier=(55.4623,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.447438,'amu*angstrom^2'), symmetry=1, barrier=(55.4623,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.447622,'amu*angstrom^2'), symmetry=1, barrier=(55.4634,'kJ/mol'), semiclassical=False),
    ],
    spinMultiplicity = 5,
    opticalIsomers = 1,
    molecularWeight = (94.1112,'amu'),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[0.381086,0.0671605,-6.36613e-05,3.14372e-08,-5.99271e-12,60162.3,28.6519], Tmin=(100,'K'), Tmax=(1428.54,'K')), NASAPolynomial(coeffs=[15.951,0.0170751,-4.25709e-06,5.35022e-10,-2.82821e-14,56376,-49.6995], Tmin=(1428.54,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(499.047,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(291.007,'J/(mol*K)'), comment="""Thermo group additivity estimation: group(O2s-CsH) + group(Cs-(Cds-Cds)OsHH) + group(Cs-(Cds-Cds)HHH) + group(Cds-CdsCsCs) + group(Cds-Cds(Cds-Cds)H) + group(Cds-Cds(Cds-Cds)H) + group(Cds-CdsHH) + radical(C=CC=CCJ) + radical(C=CCJO) + radical(C=CJC=C) + radical(C=CJC=C)"""),
)

species(
    label = '[CH][C]=C[C](C)C=O(21886)',
    structure = SMILES('[CH][C]=C[C](C)C=O'),
    E0 = (562.204,'kJ/mol'),
    modes = [
        HarmonicOscillator(frequencies=([1685,370,2750,2800,2850,1350,1500,750,1050,1375,1000,2782.5,750,1395,475,1775,1000,3010,987.5,1337.5,450,1655,360,370,350,528.116,528.158,528.386,528.436],'cm^-1')),
        HinderedRotor(inertia=(0.26683,'amu*angstrom^2'), symmetry=1, barrier=(52.7925,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.266692,'amu*angstrom^2'), symmetry=1, barrier=(52.7946,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.266691,'amu*angstrom^2'), symmetry=1, barrier=(52.7951,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.266891,'amu*angstrom^2'), symmetry=1, barrier=(52.7978,'kJ/mol'), semiclassical=False),
    ],
    spinMultiplicity = 5,
    opticalIsomers = 1,
    molecularWeight = (94.1112,'amu'),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[1.3626,0.0529383,-3.02099e-05,8.21911e-09,-8.98961e-13,67716.6,25.8835], Tmin=(100,'K'), Tmax=(2041.9,'K')), NASAPolynomial(coeffs=[14.9368,0.0263467,-1.06753e-05,1.84115e-09,-1.18069e-13,62173.2,-49.2947], Tmin=(2041.9,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(562.204,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(291.007,'J/(mol*K)'), comment="""Thermo group additivity estimation: group(Cs-(Cds-O2d)CsCsH) + group(Cs-CsHHH) + group(Cs-(Cds-Cds)HHH) + group(Cds-CdsCsH) + group(Cds-CdsCsH) + group(Cds-OdCsH) + radical(C=CCJ(C)C=O) + radical(AllylJ2_triplet) + radical(Cds_S)"""),
)

species(
    label = '[CH]C=CC([CH2])=C[O](21890)',
    structure = SMILES('[CH]C=CC([CH2])=C[O]'),
    E0 = (475.159,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
    molecularWeight = (94.1112,'amu'),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[0.573656,0.058631,-5.50884e-06,-4.6655e-08,2.55859e-11,57287.1,25.359], Tmin=(100,'K'), Tmax=(934.079,'K')), NASAPolynomial(coeffs=[20.5999,0.0152006,-3.73809e-06,5.94131e-10,-4.42022e-14,51699.4,-79.7747], Tmin=(934.079,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(475.159,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(295.164,'J/(mol*K)'), comment="""Thermo group additivity estimation: group(O2s-(Cds-Cd)H) + group(Cs-(Cds-Cds)HHH) + group(Cs-(Cds-Cds)HHH) + group(Cds-Cds(Cds-Cds)Cs) + group(Cds-CdsCsH) + group(Cds-Cds(Cds-Cds)H) + group(Cds-CdsOsH) + radical(AllylJ2_triplet) + radical(Allyl_P) + radical(C=COJ)"""),
)

species(
    label = '[CH2][C]=CC([CH2])=C[O](20121)',
    structure = SMILES('[CH2][C]=CC([CH2])=C[O]'),
    E0 = (460.372,'kJ/mol'),
    modes = [
        HarmonicOscillator(frequencies=([3000,3033.33,3066.67,3100,415,465,780,850,1435,1475,900,1100,1685,370,350,440,435,1725,2995,3025,975,1000,1300,1375,400,500,1630,1680,180,180],'cm^-1')),
        HinderedRotor(inertia=(1.67908,'amu*angstrom^2'), symmetry=1, barrier=(38.6054,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(1.67688,'amu*angstrom^2'), symmetry=1, barrier=(38.5547,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(1.67364,'amu*angstrom^2'), symmetry=1, barrier=(38.4803,'kJ/mol'), semiclassical=False),
    ],
    spinMultiplicity = 5,
    opticalIsomers = 1,
    molecularWeight = (94.1112,'amu'),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[0.81984,0.0544109,-7.43188e-06,-4.27957e-08,2.44201e-11,55498.8,25.3986], Tmin=(100,'K'), Tmax=(923.115,'K')), NASAPolynomial(coeffs=[20.2652,0.0104351,-1.43316e-06,1.46181e-10,-1.23669e-14,50192.4,-76.1553], Tmin=(923.115,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(460.372,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(295.164,'J/(mol*K)'), comment="""Thermo group additivity estimation: group(O2s-(Cds-Cd)H) + group(Cs-(Cds-Cds)HHH) + group(Cs-(Cds-Cds)HHH) + group(Cds-Cds(Cds-Cds)Cs) + group(Cds-CdsCsH) + group(Cds-Cds(Cds-Cds)H) + group(Cds-CdsOsH) + radical(Cds_S) + radical(C=CC=CCJ) + radical(C=COJ) + radical(Allyl_P)"""),
)

species(
    label = '[CH][C]=CC([CH2])=C[O](21892)',
    structure = SMILES('[CH][C]=CC([CH2])=C[O]'),
    E0 = (713,'kJ/mol'),
    modes = [
        HarmonicOscillator(frequencies=([3000,3100,440,815,1455,1000,1685,370,350,440,435,1725,2995,3025,975,1000,1300,1375,400,500,1630,1680,180,180,180,180,180],'cm^-1')),
        HinderedRotor(inertia=(2.06678,'amu*angstrom^2'), symmetry=1, barrier=(47.5194,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(2.08406,'amu*angstrom^2'), symmetry=1, barrier=(47.9166,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(2.06771,'amu*angstrom^2'), symmetry=1, barrier=(47.5407,'kJ/mol'), semiclassical=False),
    ],
    spinMultiplicity = 6,
    opticalIsomers = 1,
    molecularWeight = (93.1033,'amu'),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[0.544928,0.0631841,-3.01585e-05,-1.75119e-08,1.46432e-11,85890.3,25.9357], Tmin=(100,'K'), Tmax=(935.108,'K')), NASAPolynomial(coeffs=[19.3913,0.0147326,-4.0339e-06,6.47382e-10,-4.58793e-14,80959.3,-71.2421], Tmin=(935.108,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(713,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(270.22,'J/(mol*K)'), comment="""Thermo group additivity estimation: group(O2s-(Cds-Cd)H) + group(Cs-(Cds-Cds)HHH) + group(Cs-(Cds-Cds)HHH) + group(Cds-Cds(Cds-Cds)Cs) + group(Cds-CdsCsH) + group(Cds-Cds(Cds-Cds)H) + group(Cds-CdsOsH) + radical(Allyl_P) + radical(AllylJ2_triplet) + radical(Cds_S) + radical(C=COJ)"""),
)

species(
    label = '[CH]=[C][C]=C([CH2])[CH]O(29835)',
    structure = SMILES('[CH]=[C][C]=C([CH2])[CH]O'),
    E0 = (746.143,'kJ/mol'),
    modes = [
        HarmonicOscillator(frequencies=([3120,650,792.5,1650,3615,1277.5,1000,3025,407.5,1350,352.5,350,440,435,1725,3000,3100,440,815,1455,1000,1670,1700,300,440,361.081],'cm^-1')),
        HinderedRotor(inertia=(0.594504,'amu*angstrom^2'), symmetry=1, barrier=(55.0042,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.59451,'amu*angstrom^2'), symmetry=1, barrier=(55.0041,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.594524,'amu*angstrom^2'), symmetry=1, barrier=(55.0041,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.594513,'amu*angstrom^2'), symmetry=1, barrier=(55.004,'kJ/mol'), semiclassical=False),
    ],
    spinMultiplicity = 6,
    opticalIsomers = 1,
    molecularWeight = (93.1033,'amu'),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[0.492042,0.0685944,-7.34145e-05,4.04068e-08,-8.54353e-12,89873.9,28.6509], Tmin=(100,'K'), Tmax=(1273.82,'K')), NASAPolynomial(coeffs=[15.8097,0.014904,-3.6076e-06,4.27305e-10,-2.09677e-14,86425.1,-47.1749], Tmin=(1273.82,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(746.143,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(266.063,'J/(mol*K)'), comment="""Thermo group additivity estimation: group(O2s-CsH) + group(Cs-(Cds-Cds)OsHH) + group(Cs-(Cds-Cds)HHH) + group(Cds-CdsCsCs) + group(Cds-Cds(Cds-Cds)H) + group(Cds-Cds(Cds-Cds)H) + group(Cds-CdsHH) + radical(Cds_P) + radical(C=CC=CCJ) + radical(C=CJC=C) + radical(C=CJC=C) + radical(C=CCJO)"""),
)

species(
    label = '[CH]=[C]C1[C]([CH2])C1O(29836)',
    structure = SMILES('[CH]=[C]C1[C]([CH2])C1O'),
    E0 = (751.624,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
    molecularWeight = (94.1112,'amu'),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[1.40367,0.0447887,-5.90457e-07,-3.78717e-08,1.99295e-11,90504.3,27.4331], Tmin=(100,'K'), Tmax=(932.449,'K')), NASAPolynomial(coeffs=[15.4612,0.015757,-4.19408e-06,6.71621e-10,-4.75737e-14,86523.3,-46.6928], Tmin=(932.449,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(751.624,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(295.164,'J/(mol*K)'), comment="""Thermo group additivity estimation: group(O2s-CsH) + group(Cs-CsCsCsH) + group(Cs-(Cds-Cds)CsCsH) + group(Cs-CsCsOsH) + group(Cs-CsHHH) + group(Cds-CdsCsH) + group(Cds-CdsHH) + ring(Cyclopropane) + radical(Cds_P) + radical(Isobutyl) + radical(Cds_S) + radical(CCJ(C)CO)"""),
)

species(
    label = '[CH]=[C]C1C[C]1[CH]O(29837)',
    structure = SMILES('[CH]=[C]C1C[C]1[CH]O'),
    E0 = (739.483,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
    molecularWeight = (94.1112,'amu'),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[1.13184,0.0538741,-2.97569e-05,-3.9341e-09,6.49441e-12,89050.7,27.2141], Tmin=(100,'K'), Tmax=(984.813,'K')), NASAPolynomial(coeffs=[14.802,0.0185142,-6.61141e-06,1.18826e-09,-8.3752e-14,85380.4,-43.4917], Tmin=(984.813,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(739.483,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(295.164,'J/(mol*K)'), comment="""Thermo group additivity estimation: group(O2s-CsH) + group(Cs-CsCsCsH) + group(Cs-(Cds-Cds)CsCsH) + group(Cs-CsCsHH) + group(Cs-CsOsHH) + group(Cds-CdsCsH) + group(Cds-CdsHH) + ring(Cyclopropane) + radical(CCsJOH) + radical(Cds_S) + radical(Cds_P) + radical(CCJ(C)CO)"""),
)

species(
    label = '[CH]=C1[CH]C1([CH2])[CH]O(29797)',
    structure = SMILES('[CH]=C1[CH]C1([CH2])[CH]O'),
    E0 = (736.164,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
    molecularWeight = (94.1112,'amu'),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[0.64943,0.0651478,-5.56344e-05,2.00201e-08,-1.63388e-12,88668.1,24.9528], Tmin=(100,'K'), Tmax=(1049.99,'K')), NASAPolynomial(coeffs=[16.76,0.0173261,-6.67784e-06,1.22892e-09,-8.65554e-14,84537.8,-57.1149], Tmin=(1049.99,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(736.164,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(295.164,'J/(mol*K)'), comment="""Thermo group additivity estimation: group(O2s-CsH) + group(Cs-(Cds-Cds)CsCsCs) + group(Cs-(Cds-Cds)CsHH) + group(Cs-CsOsHH) + group(Cs-CsHHH) + group(Cds-CdsCsCs) + group(Cds-CdsHH) + ring(Methylene_cyclopropane) + radical(Neopentyl) + radical(Allyl_S) + radical(Cds_P) + radical(CCsJOH)"""),
)

species(
    label = '[CH2]C1([CH]O)[CH][C]=C1(29838)',
    structure = SMILES('[CH2]C1([CH]O)[CH][C]=C1'),
    E0 = (722.507,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
    molecularWeight = (94.1112,'amu'),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[1.05798,0.0550999,-3.23011e-05,-4.22903e-11,4.40597e-12,87011.7,23.9884], Tmin=(100,'K'), Tmax=(1038.66,'K')), NASAPolynomial(coeffs=[15.1059,0.0195804,-7.83845e-06,1.47943e-09,-1.05827e-13,83091.3,-49.1421], Tmin=(1038.66,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(722.507,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(295.164,'J/(mol*K)'), comment="""Thermo group additivity estimation: group(O2s-CsH) + group(Cs-(Cds-Cds)CsCsCs) + group(Cs-(Cds-Cds)CsHH) + group(Cs-CsOsHH) + group(Cs-CsHHH) + group(Cds-CdsCsH) + group(Cds-CdsCsH) + ring(Cyclobutene) + radical(Neopentyl) + radical(cyclobutene-vinyl) + radical(cyclobutene-allyl) + radical(CCsJOH)"""),
)

species(
    label = '[CH]=C=C=C([CH2])CO(29839)',
    structure = SMILES('[CH]=C=C=C([CH2])CO'),
    E0 = (401.066,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
    molecularWeight = (94.1112,'amu'),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[0.802665,0.0706983,-8.21366e-05,5.19693e-08,-1.31589e-11,48351.9,25.7831], Tmin=(100,'K'), Tmax=(964.028,'K')), NASAPolynomial(coeffs=[12.0735,0.0239338,-9.37362e-06,1.65161e-09,-1.10307e-13,46178.8,-28.1794], Tmin=(964.028,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(401.066,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(295.164,'J/(mol*K)'), comment="""Thermo group additivity estimation: group(O2s-CsH) + group(Cs-(Cds-Cds)OsHH) + group(Cs-(Cds-Cds)HHH) + group(Cds-CdsCsCs) + group(Cds-CdsHH) + group(Cdd-CdsCds) + group(Cdd-CdsCds) + radical(Allyl_P) + radical(C=C=CJ)"""),
)

species(
    label = '[CH]=C=C=C(C)[CH]O(29840)',
    structure = SMILES('[CH]=C=C=C(C)[CH]O'),
    E0 = (366.862,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
    molecularWeight = (94.1112,'amu'),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[0.682333,0.06604,-6.4464e-05,3.26799e-08,-6.52918e-12,44248.6,26.5288], Tmin=(100,'K'), Tmax=(1223.11,'K')), NASAPolynomial(coeffs=[15.1365,0.0187693,-6.49165e-06,1.08135e-09,-7.04752e-14,40712.9,-46.1153], Tmin=(1223.11,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(366.862,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(295.164,'J/(mol*K)'), comment="""Thermo group additivity estimation: group(O2s-CsH) + group(Cs-(Cds-Cds)OsHH) + group(Cs-(Cds-Cds)HHH) + group(Cds-CdsCsCs) + group(Cds-CdsHH) + group(Cdd-CdsCds) + group(Cdd-CdsCds) + radical(C=CCJO) + radical(C=C=CJ)"""),
)

species(
    label = '[CH]=C=C[C](C)C=O(21909)',
    structure = SMILES('[CH]=C=C[C](C)C=O'),
    E0 = (284.759,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
    molecularWeight = (94.1112,'amu'),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[1.41042,0.0498113,-2.69124e-05,1.25829e-09,2.70317e-12,34348,23.7887], Tmin=(100,'K'), Tmax=(1053.81,'K')), NASAPolynomial(coeffs=[11.3152,0.0249149,-9.55141e-06,1.7111e-09,-1.17226e-13,31555.3,-27.8608], Tmin=(1053.81,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(284.759,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(295.164,'J/(mol*K)'), comment="""Thermo group additivity estimation: group(Cs-(Cds-O2d)CsCsH) + group(Cs-CsHHH) + group(Cds-CdsCsH) + group(Cds-OdCsH) + group(Cds-CdsHH) + group(Cdd-CdsCds) + radical(C=C=CJ) + radical(C=CCJ(C)C=O)"""),
)

species(
    label = '[CH]C=CC(=C)C=O(21897)',
    structure = SMILES('[CH]C=CC(=C)C=O'),
    E0 = (314.627,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
    molecularWeight = (94.1112,'amu'),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[0.683965,0.0624036,-3.73302e-05,2.92692e-09,3.24226e-12,37969.2,24.0548], Tmin=(100,'K'), Tmax=(1083.23,'K')), NASAPolynomial(coeffs=[15.865,0.0246589,-1.04231e-05,1.97444e-09,-1.39929e-13,33605.8,-55.3583], Tmin=(1083.23,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(314.627,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(295.164,'J/(mol*K)'), comment="""Thermo group additivity estimation: group(Cs-(Cds-Cds)HHH) + group(Cds-Cds(Cds-O2d)(Cds-Cds)) + group(Cds-CdsCsH) + group(Cds-Cds(Cds-Cds)H) + group(Cds-CdsHH) + group(Cds-O2d(Cds-Cds)H) + radical(AllylJ2_triplet)"""),
)

species(
    label = 'C#C[CH][C]1CC1O(29841)',
    structure = SMILES('C#C[CH][C]1CC1O'),
    E0 = (381.89,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
    molecularWeight = (94.1112,'amu'),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[1.50545,0.0415049,1.15688e-05,-5.32295e-08,2.65842e-11,46033,23.03], Tmin=(100,'K'), Tmax=(907.9,'K')), NASAPolynomial(coeffs=[15.6689,0.0146138,-2.67058e-06,3.05987e-10,-1.98535e-14,41997.7,-51.9917], Tmin=(907.9,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(381.89,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(295.164,'J/(mol*K)'), comment="""Thermo group additivity estimation: group(O2s-CsH) + group(Cs-CsCsCsH) + group(Cs-CsCsOsH) + group(Cs-CsCsHH) + group(Cs-CtCsHH) + group(Ct-CtCs) + group(Ct-CtH) + ring(Cyclopropane) + radical(CCJ(C)CO) + radical(Sec_Propargyl)"""),
)

species(
    label = '[CH]=C1[CH]C(=C)C1O(29795)',
    structure = SMILES('[CH]=C1[CH]C(=C)C1O'),
    E0 = (367.124,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
    molecularWeight = (94.1112,'amu'),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[1.89318,0.0223674,7.7443e-05,-1.24941e-07,5.16925e-11,44252.7,23.7356], Tmin=(100,'K'), Tmax=(939.361,'K')), NASAPolynomial(coeffs=[18.9084,0.0108998,-1.63039e-06,2.91829e-10,-3.08236e-14,38365.3,-71.6105], Tmin=(939.361,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(367.124,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(303.478,'J/(mol*K)'), comment="""Thermo group additivity estimation: group(O2s-CsH) + group(Cs-(Cds-Cds)(Cds-Cds)OsH) + group(Cs-(Cds-Cds)(Cds-Cds)HH) + group(Cds-CdsCsCs) + group(Cds-CdsCsCs) + group(Cds-CdsHH) + group(Cds-CdsHH) + ring(Cyclobutane) + radical(C=CCJC=C) + radical(Cds_P)"""),
)

species(
    label = '[CH]=C1[CH]C(=CO)C1(29790)',
    structure = SMILES('[CH]=C1[CH]C(=CO)C1'),
    E0 = (336.165,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
    molecularWeight = (94.1112,'amu'),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[1.5366,0.0213298,0.000106852,-1.76472e-07,7.59168e-11,40550.6,22.4922], Tmin=(100,'K'), Tmax=(914.296,'K')), NASAPolynomial(coeffs=[26.7533,-0.00345981,7.19759e-06,-1.48919e-09,9.28454e-14,32364.5,-116.455], Tmin=(914.296,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(336.165,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(303.478,'J/(mol*K)'), comment="""Thermo group additivity estimation: group(O2s-(Cds-Cd)H) + group(Cs-(Cds-Cds)(Cds-Cds)HH) + group(Cs-(Cds-Cds)(Cds-Cds)HH) + group(Cds-CdsCsCs) + group(Cds-CdsCsCs) + group(Cds-CdsOsH) + group(Cds-CdsHH) + ring(Cyclobutane) + radical(C=CCJC=C) + radical(Cds_P)"""),
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
    label = '[CH]C([CH2])=C[C]=[CH](19887)',
    structure = SMILES('[CH]C([CH2])=C[C]=[CH]'),
    E0 = (953.672,'kJ/mol'),
    modes = [
        HarmonicOscillator(frequencies=([3000,3100,440,815,1455,1000,1685,370,350,440,435,1725,3120,650,792.5,1650,3010,987.5,1337.5,450,1655,180,180,180],'cm^-1')),
        HinderedRotor(inertia=(2.198,'amu*angstrom^2'), symmetry=1, barrier=(50.5364,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(2.19802,'amu*angstrom^2'), symmetry=1, barrier=(50.5367,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(2.1981,'amu*angstrom^2'), symmetry=1, barrier=(50.5387,'kJ/mol'), semiclassical=False),
    ],
    spinMultiplicity = 6,
    opticalIsomers = 1,
    molecularWeight = (77.1039,'amu'),
    collisionModel = TransportData(shapeIndex=2, epsilon=(3312.87,'J/mol'), sigma=(5.74688,'angstroms'), dipoleMoment=(0,'C*m'), polarizability=(0,'angstroms^3'), rotrelaxcollnum=0, comment="""Epsilon & sigma estimated with Tc=517.46 K, Pc=39.61 bar (from Joback method)"""),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[1.33019,0.0534685,-3.53423e-05,5.78722e-09,2.75888e-12,114801,22.8075], Tmin=(100,'K'), Tmax=(930.983,'K')), NASAPolynomial(coeffs=[11.325,0.0235591,-8.15178e-06,1.354e-09,-8.87473e-14,112375,-27.7303], Tmin=(930.983,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(953.672,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(245.277,'J/(mol*K)'), comment="""Thermo group additivity estimation: group(Cs-(Cds-Cds)HHH) + group(Cs-(Cds-Cds)HHH) + group(Cds-CdsCsCs) + group(Cds-Cds(Cds-Cds)H) + group(Cds-Cds(Cds-Cds)H) + group(Cds-CdsHH) + radical(C=CC=CCJ) + radical(C=CJC=C) + radical(Cds_P) + radical(AllylJ2_triplet)"""),
)

species(
    label = '[CH]=[C]C=[C][CH2](19274)',
    structure = SMILES('[CH]=[C]C=[C][CH2]'),
    E0 = (859.885,'kJ/mol'),
    modes = [
        HarmonicOscillator(frequencies=([3120,650,792.5,1650,1670,1700,300,440,3010,987.5,1337.5,450,1655,3000,3100,440,815,1455,1000],'cm^-1')),
        HinderedRotor(inertia=(2.0544,'amu*angstrom^2'), symmetry=1, barrier=(47.2346,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(2.05464,'amu*angstrom^2'), symmetry=1, barrier=(47.2402,'kJ/mol'), semiclassical=False),
    ],
    spinMultiplicity = 5,
    opticalIsomers = 1,
    molecularWeight = (64.0853,'amu'),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[1.93557,0.0405981,-3.87335e-05,2.02214e-08,-4.17345e-12,103499,19.2405], Tmin=(100,'K'), Tmax=(1265.13,'K')), NASAPolynomial(coeffs=[9.93801,0.0132358,-3.8481e-06,5.50841e-10,-3.19634e-14,101639,-20.5969], Tmin=(1265.13,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(859.885,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(199.547,'J/(mol*K)'), comment="""Thermo group additivity estimation: group(Cs-(Cds-Cds)HHH) + group(Cds-CdsCsH) + group(Cds-Cds(Cds-Cds)H) + group(Cds-Cds(Cds-Cds)H) + group(Cds-CdsHH) + radical(Cds_S) + radical(C=CC=CCJ) + radical(C=CJC=C) + radical(Cds_P)"""),
)

species(
    label = '[CH][C]=C[C]=CO(29842)',
    structure = SMILES('[CH][C]=C[C]=CO'),
    E0 = (656.625,'kJ/mol'),
    modes = [
        HarmonicOscillator(frequencies=([1670,1700,300,440,3615,1277.5,1000,2995,3025,975,1000,1300,1375,400,500,1630,1680,180,180,180,180],'cm^-1')),
        HinderedRotor(inertia=(2.01369,'amu*angstrom^2'), symmetry=1, barrier=(46.2988,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(2.01047,'amu*angstrom^2'), symmetry=1, barrier=(46.2246,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(2.01467,'amu*angstrom^2'), symmetry=1, barrier=(46.3212,'kJ/mol'), semiclassical=False),
    ],
    spinMultiplicity = 5,
    opticalIsomers = 1,
    molecularWeight = (80.0847,'amu'),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[0.120588,0.0683853,-7.29522e-05,3.76676e-08,-7.20697e-12,79127.9,24.763], Tmin=(100,'K'), Tmax=(1505,'K')), NASAPolynomial(coeffs=[18.1039,0.00788404,1.08304e-08,-2.61966e-10,2.53755e-14,75153.8,-64.5672], Tmin=(1505,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(656.625,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(220.334,'J/(mol*K)'), comment="""Thermo group additivity estimation: group(O2s-(Cds-Cd)H) + group(Cs-(Cds-Cds)HHH) + group(Cds-CdsCsH) + group(Cds-Cds(Cds-Cds)H) + group(Cds-Cds(Cds-Cds)H) + group(Cds-CdsOsH) + radical(AllylJ2_triplet) + radical(Cds_S) + radical(C=CJC=C)"""),
)

species(
    label = '[CH]C([CH2])=CO(19098)',
    structure = SMILES('[CH]C([CH2])=CO'),
    E0 = (280.46,'kJ/mol'),
    modes = [
        HarmonicOscillator(frequencies=([3000,3100,440,815,1455,1000,3010,987.5,1337.5,450,1655,3615,1277.5,1000,350,440,435,1725,180,180,180],'cm^-1')),
        HinderedRotor(inertia=(2.07603,'amu*angstrom^2'), symmetry=1, barrier=(47.7321,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(2.07698,'amu*angstrom^2'), symmetry=1, barrier=(47.7538,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(2.07646,'amu*angstrom^2'), symmetry=1, barrier=(47.7419,'kJ/mol'), semiclassical=False),
    ],
    spinMultiplicity = 4,
    opticalIsomers = 1,
    molecularWeight = (69.0819,'amu'),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[1.45652,0.0411178,9.02269e-06,-5.40178e-08,2.77613e-11,33836.9,18.3017], Tmin=(100,'K'), Tmax=(910.922,'K')), NASAPolynomial(coeffs=[18.0563,0.00742605,-4.82671e-08,-1.37232e-10,8.52061e-15,29186.3,-69.1609], Tmin=(910.922,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(280.46,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(220.334,'J/(mol*K)'), comment="""Thermo group additivity estimation: group(O2s-(Cds-Cd)H) + group(Cs-(Cds-Cds)HHH) + group(Cs-(Cds-Cds)HHH) + group(Cds-CdsCsCs) + group(Cds-CdsOsH) + radical(AllylJ2_triplet) + radical(Allyl_P)"""),
)

species(
    label = '[C]=[CH](18830)',
    structure = SMILES('[C]=[CH]'),
    E0 = (847.348,'kJ/mol'),
    modes = [
        HarmonicOscillator(frequencies=([562.459,1392.74,3112.83],'cm^-1')),
    ],
    spinMultiplicity = 4,
    opticalIsomers = 1,
    molecularWeight = (25.0293,'amu'),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[3.86574,0.0023643,1.31489e-06,-2.64796e-09,1.00932e-12,101918,6.62295], Tmin=(100,'K'), Tmax=(1016.96,'K')), NASAPolynomial(coeffs=[4.17915,0.00251347,-9.43473e-07,1.6872e-10,-1.15831e-14,101782,4.75427], Tmin=(1016.96,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(847.348,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(58.2013,'J/(mol*K)'), comment="""Thermo group additivity estimation: group(Cds-CdsHH) + group(Cds-CdsHH) + radical(Cds_P) + radical(CdCdJ2_triplet)"""),
)

species(
    label = '[CH][C]=CC([CH2])=[C]O(29843)',
    structure = SMILES('[CH][C]=CC([CH2])=[C]O'),
    E0 = (811.282,'kJ/mol'),
    modes = [
        HarmonicOscillator(frequencies=([3000,3100,440,815,1455,1000,1670,1700,300,440,3010,987.5,1337.5,450,1655,3615,1277.5,1000,350,440,435,1725,180,180,180,180],'cm^-1')),
        HinderedRotor(inertia=(2.08472,'amu*angstrom^2'), symmetry=1, barrier=(47.9319,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(2.08518,'amu*angstrom^2'), symmetry=1, barrier=(47.9425,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(2.08638,'amu*angstrom^2'), symmetry=1, barrier=(47.9699,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(2.08599,'amu*angstrom^2'), symmetry=1, barrier=(47.9611,'kJ/mol'), semiclassical=False),
    ],
    spinMultiplicity = 6,
    opticalIsomers = 1,
    molecularWeight = (93.1033,'amu'),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[-0.162921,0.0778568,-8.11274e-05,4.17473e-08,-8.1754e-12,97736.3,30.1956], Tmin=(100,'K'), Tmax=(1384.11,'K')), NASAPolynomial(coeffs=[19.8035,0.012947,-2.97104e-06,3.40204e-10,-1.67776e-14,92899.6,-70.1268], Tmin=(1384.11,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(811.282,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(266.063,'J/(mol*K)'), comment="""Thermo group additivity estimation: group(O2s-(Cds-Cd)H) + group(Cs-(Cds-Cds)HHH) + group(Cs-(Cds-Cds)HHH) + group(Cds-Cds(Cds-Cds)Cs) + group(Cds-CdsCsH) + group(Cds-Cds(Cds-Cds)H) + group(Cds-CdsOsH) + radical(C=CJO) + radical(Allyl_P) + radical(AllylJ2_triplet) + radical(Cds_S)"""),
)

species(
    label = '[CH][C]=CC([CH])=CO(24780)',
    structure = SMILES('[CH][C]=CC([CH])=CO'),
    E0 = (790.723,'kJ/mol'),
    modes = [
        HarmonicOscillator(frequencies=([1685,370,2995,3025,975,1000,1300,1375,400,500,1630,1680,3615,1277.5,1000,350,440,435,1725,200,800,960,1120,1280,1440,1600],'cm^-1')),
        HinderedRotor(inertia=(0.156089,'amu*angstrom^2'), symmetry=1, barrier=(3.5888,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.156089,'amu*angstrom^2'), symmetry=1, barrier=(3.5888,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.156089,'amu*angstrom^2'), symmetry=1, barrier=(3.5888,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.156089,'amu*angstrom^2'), symmetry=1, barrier=(3.5888,'kJ/mol'), semiclassical=False),
    ],
    spinMultiplicity = 6,
    opticalIsomers = 1,
    molecularWeight = (93.1033,'amu'),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[-0.84554,0.0838804,-8.02754e-05,3.78343e-08,-6.73121e-12,95295.9,30.3604], Tmin=(100,'K'), Tmax=(1583.5,'K')), NASAPolynomial(coeffs=[21.4952,0.0148752,-3.00051e-06,2.87324e-10,-1.18326e-14,89796.7,-82.7129], Tmin=(1583.5,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(790.723,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(266.063,'J/(mol*K)'), comment="""Thermo group additivity estimation: group(O2s-(Cds-Cd)H) + group(Cs-(Cds-Cds)HHH) + group(Cs-(Cds-Cds)HHH) + group(Cds-Cds(Cds-Cds)Cs) + group(Cds-CdsCsH) + group(Cds-Cds(Cds-Cds)H) + group(Cds-CdsOsH) + radical(AllylJ2_triplet) + radical(AllylJ2_triplet) + radical(Cds_S)"""),
)

species(
    label = '[C]=[C]C=C([CH2])[CH]O(29844)',
    structure = SMILES('[C]=[C]C=C([CH2])[CH]O'),
    E0 = (858.153,'kJ/mol'),
    modes = [
        HarmonicOscillator(frequencies=([3615,1277.5,1000,3025,407.5,1350,352.5,350,440,435,1725,3000,3100,440,815,1455,1000,1685,370,3010,987.5,1337.5,450,1655,469.621,469.623],'cm^-1')),
        HinderedRotor(inertia=(0.539172,'amu*angstrom^2'), symmetry=1, barrier=(84.3829,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.539173,'amu*angstrom^2'), symmetry=1, barrier=(84.3829,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.539174,'amu*angstrom^2'), symmetry=1, barrier=(84.3828,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.539176,'amu*angstrom^2'), symmetry=1, barrier=(84.3829,'kJ/mol'), semiclassical=False),
    ],
    spinMultiplicity = 6,
    opticalIsomers = 1,
    molecularWeight = (93.1033,'amu'),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[0.589218,0.0650237,-6.48237e-05,3.28724e-08,-6.44707e-12,103343,28.6862], Tmin=(100,'K'), Tmax=(1336.45,'K')), NASAPolynomial(coeffs=[16.2847,0.0144848,-4.10179e-06,5.87805e-10,-3.4739e-14,99466.2,-50.3972], Tmin=(1336.45,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(858.153,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(266.063,'J/(mol*K)'), comment="""Thermo group additivity estimation: group(O2s-CsH) + group(Cs-(Cds-Cds)OsHH) + group(Cs-(Cds-Cds)HHH) + group(Cds-CdsCsCs) + group(Cds-Cds(Cds-Cds)H) + group(Cds-Cds(Cds-Cds)H) + group(Cds-CdsHH) + radical(C=CJC=C) + radical(C=CCJO) + radical(C=CC=CCJ) + radical(CdCdJ2_triplet)"""),
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
    label = '[CH][C]=CC(C)=[C]O(29845)',
    structure = SMILES('[CH][C]=CC(C)=[C]O'),
    E0 = (659.783,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
    molecularWeight = (94.1112,'amu'),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[0.322578,0.0746707,-7.00233e-05,2.99336e-08,-3.59069e-12,79491.6,28.1709], Tmin=(100,'K'), Tmax=(949.702,'K')), NASAPolynomial(coeffs=[16.38,0.0209378,-7.10738e-06,1.17804e-09,-7.75581e-14,75814.8,-51.7682], Tmin=(949.702,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(659.783,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(291.007,'J/(mol*K)'), comment="""Thermo group additivity estimation: group(O2s-(Cds-Cd)H) + group(Cs-(Cds-Cds)HHH) + group(Cs-(Cds-Cds)HHH) + group(Cds-Cds(Cds-Cds)Cs) + group(Cds-CdsCsH) + group(Cds-Cds(Cds-Cds)H) + group(Cds-CdsOsH) + radical(C=CJO) + radical(AllylJ2_triplet) + radical(Cds_S)"""),
)

species(
    label = '[CH]C=CC([CH2])=[C]O(29846)',
    structure = SMILES('[CH]C=CC([CH2])=[C]O'),
    E0 = (573.44,'kJ/mol'),
    modes = [
        HarmonicOscillator(frequencies=([3000,3100,440,815,1455,1000,1685,370,2995,3025,975,1000,1300,1375,400,500,1630,1680,3615,1277.5,1000,350,440,435,1725,180,180,180,180],'cm^-1')),
        HinderedRotor(inertia=(2.07226,'amu*angstrom^2'), symmetry=1, barrier=(47.6454,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(2.07229,'amu*angstrom^2'), symmetry=1, barrier=(47.6459,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(2.07222,'amu*angstrom^2'), symmetry=1, barrier=(47.6445,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(2.07228,'amu*angstrom^2'), symmetry=1, barrier=(47.6458,'kJ/mol'), semiclassical=False),
    ],
    spinMultiplicity = 5,
    opticalIsomers = 1,
    molecularWeight = (94.1112,'amu'),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[0.37397,0.0673299,-3.56596e-05,-1.43206e-08,1.42393e-11,69111,27.7947], Tmin=(100,'K'), Tmax=(923.16,'K')), NASAPolynomial(coeffs=[19.7624,0.0155569,-3.91461e-06,5.80114e-10,-3.94261e-14,64157.6,-71.6328], Tmin=(923.16,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(573.44,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(291.007,'J/(mol*K)'), comment="""Thermo group additivity estimation: group(O2s-(Cds-Cd)H) + group(Cs-(Cds-Cds)HHH) + group(Cs-(Cds-Cds)HHH) + group(Cds-Cds(Cds-Cds)Cs) + group(Cds-CdsCsH) + group(Cds-Cds(Cds-Cds)H) + group(Cds-CdsOsH) + radical(C=CJO) + radical(AllylJ2_triplet) + radical(Allyl_P)"""),
)

species(
    label = '[CH][C]=C[C]1CC1O(29847)',
    structure = SMILES('[CH][C]=C[C]1CC1O'),
    E0 = (666.927,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
    molecularWeight = (94.1112,'amu'),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[1.45664,0.0426941,1.26331e-05,-4.62532e-08,2.05409e-11,80316.1,26.7014], Tmin=(100,'K'), Tmax=(978.859,'K')), NASAPolynomial(coeffs=[13.3565,0.0256058,-9.51105e-06,1.74449e-09,-1.24464e-13,76475.5,-38.1722], Tmin=(978.859,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(666.927,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(295.164,'J/(mol*K)'), comment="""Thermo group additivity estimation: group(O2s-CsH) + group(Cs-(Cds-Cds)CsCsH) + group(Cs-CsCsOsH) + group(Cs-CsCsHH) + group(Cs-(Cds-Cds)HHH) + group(Cds-CdsCsH) + group(Cds-CdsCsH) + ring(Cyclopropane) + radical(CCJ(C)CO) + radical(Cds_S) + radical(AllylJ2_triplet)"""),
)

species(
    label = '[CH]C1=C[C]([CH2])C1O(29848)',
    structure = SMILES('[CH]C1=C[C]([CH2])C1O'),
    E0 = (629.194,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
    molecularWeight = (94.1112,'amu'),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[1.49655,0.0382575,3.38389e-05,-7.48666e-08,3.26542e-11,75780.1,27.3946], Tmin=(100,'K'), Tmax=(942.691,'K')), NASAPolynomial(coeffs=[15.3906,0.0214655,-6.53077e-06,1.12765e-09,-8.16414e-14,71287.1,-48.7527], Tmin=(942.691,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(629.194,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(295.164,'J/(mol*K)'), comment="""Thermo group additivity estimation: group(O2s-CsH) + group(Cs-(Cds-Cds)CsCsH) + group(Cs-(Cds-Cds)CsOsH) + group(Cs-CsHHH) + group(Cs-(Cds-Cds)HHH) + group(Cds-CdsCsCs) + group(Cds-CdsCsH) + ring(Cyclobutene) + radical(CCJ(C)CO) + radical(AllylJ2_triplet) + radical(Isobutyl)"""),
)

species(
    label = '[CH][C]=C[C]([CH2])C[O](21898)',
    structure = SMILES('[CH][C]=C[C]([CH2])C[O]'),
    E0 = (947.761,'kJ/mol'),
    modes = [
        HarmonicOscillator(frequencies=([3000,3100,440,815,1455,1000,2750,2850,1437.5,1250,1305,750,350,3010,987.5,1337.5,450,1655,1685,370,360,370,350,200,800,1000,1200,1400,1600],'cm^-1')),
        HinderedRotor(inertia=(0.156089,'amu*angstrom^2'), symmetry=1, barrier=(3.5888,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.156089,'amu*angstrom^2'), symmetry=1, barrier=(3.5888,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.156089,'amu*angstrom^2'), symmetry=1, barrier=(3.5888,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.156089,'amu*angstrom^2'), symmetry=1, barrier=(3.5888,'kJ/mol'), semiclassical=False),
    ],
    spinMultiplicity = 7,
    opticalIsomers = 1,
    molecularWeight = (94.1112,'amu'),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[1.34147,0.0634599,-7.19899e-05,5.91174e-08,-2.108e-11,114080,31.2443], Tmin=(100,'K'), Tmax=(813.352,'K')), NASAPolynomial(coeffs=[3.88481,0.0428668,-1.91009e-05,3.54512e-09,-2.42113e-13,113934,21.1436], Tmin=(813.352,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(947.761,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(291.007,'J/(mol*K)'), comment="""Thermo group additivity estimation: group(O2s-CsH) + group(Cs-(Cds-Cds)CsCsH) + group(Cs-CsOsHH) + group(Cs-CsHHH) + group(Cs-(Cds-Cds)HHH) + group(Cds-CdsCsH) + group(Cds-CdsCsH) + radical(CCOJ) + radical(Cds_S) + radical(AllylJ2_triplet) + radical(Isobutyl) + radical(CCJ(C)CO)"""),
)

species(
    label = '[CH][C]=[C]C([CH2])[CH]O(21900)',
    structure = SMILES('[CH][C]=[C]C([CH2])[CH]O'),
    E0 = (987.621,'kJ/mol'),
    modes = [
        HarmonicOscillator(frequencies=([3000,3100,440,815,1455,1000,3025,407.5,1350,352.5,3615,1277.5,1000,1670,1700,300,440,1380,1390,370,380,2900,435,200,800,1066.67,1333.33,1600],'cm^-1')),
        HinderedRotor(inertia=(0.156089,'amu*angstrom^2'), symmetry=1, barrier=(3.5888,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.156089,'amu*angstrom^2'), symmetry=1, barrier=(3.5888,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.156089,'amu*angstrom^2'), symmetry=1, barrier=(3.5888,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.156089,'amu*angstrom^2'), symmetry=1, barrier=(3.5888,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.156089,'amu*angstrom^2'), symmetry=1, barrier=(3.5888,'kJ/mol'), semiclassical=False),
    ],
    spinMultiplicity = 7,
    opticalIsomers = 1,
    molecularWeight = (94.1112,'amu'),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[0.382101,0.0864791,-0.000128028,1.06466e-07,-3.4597e-11,118907,32.9603], Tmin=(100,'K'), Tmax=(894.37,'K')), NASAPolynomial(coeffs=[8.63256,0.0341729,-1.44624e-05,2.5528e-09,-1.66677e-13,118048,-2.47753], Tmin=(894.37,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(987.621,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(286.849,'J/(mol*K)'), comment="""Thermo group additivity estimation: group(O2s-CsH) + group(Cs-(Cds-Cds)CsCsH) + group(Cs-CsOsHH) + group(Cs-CsHHH) + group(Cs-(Cds-Cds)HHH) + group(Cds-CdsCsH) + group(Cds-CdsCsH) + radical(Isobutyl) + radical(Cds_S) + radical(AllylJ2_triplet) + radical(Cds_S) + radical(CCsJOH)"""),
)

species(
    label = '[CH][C]=CC[C]=CO(21983)',
    structure = SMILES('[CH][C]=CC[C]=CO'),
    E0 = (690.3,'kJ/mol'),
    modes = [
        HarmonicOscillator(frequencies=([1670,1700,300,440,2750,2850,1437.5,1250,1305,750,350,3615,1277.5,1000,2995,3025,975,1000,1300,1375,400,500,1630,1680,200,800,1066.67,1333.33,1600],'cm^-1')),
        HinderedRotor(inertia=(0.156089,'amu*angstrom^2'), symmetry=1, barrier=(3.5888,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.156089,'amu*angstrom^2'), symmetry=1, barrier=(3.5888,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.156089,'amu*angstrom^2'), symmetry=1, barrier=(3.5888,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.156089,'amu*angstrom^2'), symmetry=1, barrier=(3.5888,'kJ/mol'), semiclassical=False),
    ],
    spinMultiplicity = 5,
    opticalIsomers = 1,
    molecularWeight = (94.1112,'amu'),
    collisionModel = TransportData(shapeIndex=2, epsilon=(3932.86,'J/mol'), sigma=(6.48948,'angstroms'), dipoleMoment=(0,'C*m'), polarizability=(0,'angstroms^3'), rotrelaxcollnum=0, comment="""Epsilon & sigma estimated with Tc=614.30 K, Pc=32.65 bar (from Joback method)"""),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[0.475394,0.0686749,-5.1018e-05,9.47526e-09,3.69723e-12,83158.8,29.1952], Tmin=(100,'K'), Tmax=(951.975,'K')), NASAPolynomial(coeffs=[16.9802,0.0199179,-6.64032e-06,1.12049e-09,-7.6002e-14,79083.2,-54.5197], Tmin=(951.975,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(690.3,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(291.007,'J/(mol*K)'), comment="""Thermo group additivity estimation: group(O2s-(Cds-Cd)H) + group(Cs-(Cds-Cds)(Cds-Cds)HH) + group(Cs-(Cds-Cds)HHH) + group(Cds-CdsCsH) + group(Cds-CdsCsH) + group(Cds-CdsCsH) + group(Cds-CdsOsH) + radical(Cds_S) + radical(AllylJ2_triplet) + radical(Cds_S)"""),
)

species(
    label = '[CH][C]=CC([CH2])C=O(19974)',
    structure = SMILES('[CH][C]=CC([CH2])C=O'),
    E0 = (653.015,'kJ/mol'),
    modes = [
        HarmonicOscillator(frequencies=([3000,3100,440,815,1455,1000,2782.5,750,1395,475,1775,1000,1380,1390,370,380,2900,435,1685,370,3010,987.5,1337.5,450,1655,418.288,419.792,421.666,423.663],'cm^-1')),
        HinderedRotor(inertia=(0.432176,'amu*angstrom^2'), symmetry=1, barrier=(54.3801,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.425636,'amu*angstrom^2'), symmetry=1, barrier=(54.446,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.425028,'amu*angstrom^2'), symmetry=1, barrier=(54.3678,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.43359,'amu*angstrom^2'), symmetry=1, barrier=(54.3871,'kJ/mol'), semiclassical=False),
    ],
    spinMultiplicity = 5,
    opticalIsomers = 1,
    molecularWeight = (94.1112,'amu'),
    collisionModel = TransportData(shapeIndex=2, epsilon=(3682.09,'J/mol'), sigma=(6.19766,'angstroms'), dipoleMoment=(0,'C*m'), polarizability=(0,'angstroms^3'), rotrelaxcollnum=0, comment="""Epsilon & sigma estimated with Tc=575.13 K, Pc=35.1 bar (from Joback method)"""),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[0.807588,0.0774699,-0.000108726,9.4887e-08,-3.36644e-11,78647.5,28.7471], Tmin=(100,'K'), Tmax=(823.805,'K')), NASAPolynomial(coeffs=[5.65568,0.0411331,-1.92632e-05,3.63267e-09,-2.49176e-13,78282.9,8.93298], Tmin=(823.805,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(653.015,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(291.007,'J/(mol*K)'), comment="""Thermo group additivity estimation: group(Cs-(Cds-O2d)CsCsH) + group(Cs-CsHHH) + group(Cs-(Cds-Cds)HHH) + group(Cds-CdsCsH) + group(Cds-CdsCsH) + group(Cds-OdCsH) + radical(AllylJ2_triplet) + radical(Cds_S) + radical(CJC(C)C=O)"""),
)

species(
    label = '[CH]=[C]CC([CH2])=[C]O(28748)',
    structure = SMILES('[CH]=[C]CC([CH2])=[C]O'),
    E0 = (717.083,'kJ/mol'),
    modes = [
        HarmonicOscillator(frequencies=([3120,650,792.5,1650,3615,1277.5,1000,350,440,435,1725,2750,2850,1437.5,1250,1305,750,350,3000,3100,440,815,1455,1000,1670,1700,300,440,302.952],'cm^-1')),
        HinderedRotor(inertia=(0.266361,'amu*angstrom^2'), symmetry=1, barrier=(17.3873,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.265561,'amu*angstrom^2'), symmetry=1, barrier=(17.3848,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.265414,'amu*angstrom^2'), symmetry=1, barrier=(17.3864,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.267776,'amu*angstrom^2'), symmetry=1, barrier=(17.3882,'kJ/mol'), semiclassical=False),
    ],
    spinMultiplicity = 5,
    opticalIsomers = 1,
    molecularWeight = (94.1112,'amu'),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[0.229545,0.0739095,-8.03635e-05,4.33507e-08,-8.98668e-12,86388.7,31.2584], Tmin=(100,'K'), Tmax=(1231.1,'K')), NASAPolynomial(coeffs=[18.4045,0.0129135,-3.6769e-06,5.4115e-10,-3.2963e-14,82060.9,-59.6061], Tmin=(1231.1,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(717.083,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(291.007,'J/(mol*K)'), comment="""Thermo group additivity estimation: group(O2s-(Cds-Cd)H) + group(Cs-(Cds-Cds)(Cds-Cds)HH) + group(Cs-(Cds-Cds)HHH) + group(Cds-CdsCsCs) + group(Cds-CdsCsH) + group(Cds-CdsOsH) + group(Cds-CdsHH) + radical(Allyl_P) + radical(Cds_S) + radical(Cds_P) + radical(C=CJO)"""),
)

species(
    label = '[CH]=[C]CC([CH2])=C[O](19973)',
    structure = SMILES('[CH]=[C]CC([CH2])=C[O]'),
    E0 = (618.802,'kJ/mol'),
    modes = [
        HarmonicOscillator(frequencies=([3120,650,792.5,1650,350,440,435,1725,2750,2850,1437.5,1250,1305,750,350,3000,3100,440,815,1455,1000,1685,370,3010,987.5,1337.5,450,1655,273.131,273.81],'cm^-1')),
        HinderedRotor(inertia=(0.425078,'amu*angstrom^2'), symmetry=1, barrier=(22.8387,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.424101,'amu*angstrom^2'), symmetry=1, barrier=(22.8304,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.425145,'amu*angstrom^2'), symmetry=1, barrier=(22.8309,'kJ/mol'), semiclassical=False),
    ],
    spinMultiplicity = 5,
    opticalIsomers = 1,
    molecularWeight = (94.1112,'amu'),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[0.676551,0.0622972,-4.00564e-05,-2.07874e-09,7.90601e-12,74554.1,27.9356], Tmin=(100,'K'), Tmax=(957.806,'K')), NASAPolynomial(coeffs=[18.2744,0.0141737,-4.42049e-06,7.70418e-10,-5.54612e-14,70019.4,-62.2795], Tmin=(957.806,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(618.802,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(295.164,'J/(mol*K)'), comment="""Thermo group additivity estimation: group(O2s-(Cds-Cd)H) + group(Cs-(Cds-Cds)(Cds-Cds)HH) + group(Cs-(Cds-Cds)HHH) + group(Cds-CdsCsCs) + group(Cds-CdsCsH) + group(Cds-CdsOsH) + group(Cds-CdsHH) + radical(Cds_S) + radical(Cds_P) + radical(C=COJ) + radical(Allyl_P)"""),
)

species(
    label = '[CH2][C]=CC([CH2])=[C]O(28971)',
    structure = SMILES('[CH2][C]=CC([CH2])=[C]O'),
    E0 = (558.653,'kJ/mol'),
    modes = [
        HarmonicOscillator(frequencies=([3000,3033.33,3066.67,3100,415,465,780,850,1435,1475,900,1100,1670,1700,300,440,3010,987.5,1337.5,450,1655,3615,1277.5,1000,350,440,435,1725,226.963],'cm^-1')),
        HinderedRotor(inertia=(0.709261,'amu*angstrom^2'), symmetry=1, barrier=(25.967,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.705851,'amu*angstrom^2'), symmetry=1, barrier=(25.977,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.696119,'amu*angstrom^2'), symmetry=1, barrier=(26.0068,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(1.64168,'amu*angstrom^2'), symmetry=1, barrier=(61.1827,'kJ/mol'), semiclassical=False),
    ],
    spinMultiplicity = 5,
    opticalIsomers = 1,
    molecularWeight = (94.1112,'amu'),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[-0.289585,0.0736764,-7.3901e-05,3.59471e-08,-6.50941e-12,67362.7,31.11], Tmin=(100,'K'), Tmax=(1579.12,'K')), NASAPolynomial(coeffs=[19.9684,0.00951018,-7.41996e-07,-9.26983e-11,1.21963e-14,62567.1,-70.8046], Tmin=(1579.12,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(558.653,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(291.007,'J/(mol*K)'), comment="""Thermo group additivity estimation: group(O2s-(Cds-Cd)H) + group(Cs-(Cds-Cds)HHH) + group(Cs-(Cds-Cds)HHH) + group(Cds-Cds(Cds-Cds)Cs) + group(Cds-CdsCsH) + group(Cds-Cds(Cds-Cds)H) + group(Cds-CdsOsH) + radical(Cds_S) + radical(C=CJO) + radical(Allyl_P) + radical(C=CC=CCJ)"""),
)

species(
    label = '[CH2][C]1C=[C][CH]C1O(29849)',
    structure = SMILES('[CH2][C]1C=[C][CH]C1O'),
    E0 = (549.854,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
    molecularWeight = (94.1112,'amu'),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[1.5939,0.0368393,2.77874e-05,-6.69427e-08,2.95965e-11,66233.4,23.308], Tmin=(100,'K'), Tmax=(948.283,'K')), NASAPolynomial(coeffs=[15.6837,0.0175616,-5.23703e-06,9.29111e-10,-6.94191e-14,61755.8,-53.4387], Tmin=(948.283,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(549.854,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(299.321,'J/(mol*K)'), comment="""Thermo group additivity estimation: group(O2s-CsH) + group(Cs-(Cds-Cds)CsCsH) + group(Cs-CsCsOsH) + group(Cs-(Cds-Cds)CsHH) + group(Cs-CsHHH) + group(Cds-CdsCsH) + group(Cds-CdsCsH) + ring(Cyclopentene) + radical(Isobutyl) + radical(C=CCJCO) + radical(CCJ(C)CO) + radical(cyclopentene-vinyl)"""),
)

species(
    label = '[CH]=[C]C1CC1=CO(29817)',
    structure = SMILES('[CH]=[C]C1CC1=CO'),
    E0 = (534.791,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
    molecularWeight = (94.1112,'amu'),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[-0.365981,0.0741274,-7.47162e-05,3.59312e-08,-6.40894e-12,64496.6,26.2383], Tmin=(100,'K'), Tmax=(1603.21,'K')), NASAPolynomial(coeffs=[20.9908,0.0076964,-2.62231e-07,-1.43795e-10,1.40611e-14,59338.1,-81.6076], Tmin=(1603.21,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(534.791,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(299.321,'J/(mol*K)'), comment="""Thermo group additivity estimation: group(O2s-(Cds-Cd)H) + group(Cs-(Cds-Cds)(Cds-Cds)CsH) + group(Cs-(Cds-Cds)CsHH) + group(Cds-CdsCsCs) + group(Cds-CdsCsH) + group(Cds-CdsOsH) + group(Cds-CdsHH) + ring(Methylene_cyclopropane) + radical(Cds_P) + radical(Cds_S)"""),
)

species(
    label = '[CH][C]([CH]C#C)CO(29850)',
    structure = SMILES('[CH][C]([CH]C#C)CO'),
    E0 = (680.151,'kJ/mol'),
    modes = [
        HarmonicOscillator(frequencies=([360,370,350,750,770,3400,2100,3025,407.5,1350,352.5,2750,2850,1437.5,1250,1305,750,350,2175,525,3615,1277.5,1000,200,800,1066.67,1333.33,1600],'cm^-1')),
        HinderedRotor(inertia=(0.156089,'amu*angstrom^2'), symmetry=1, barrier=(3.5888,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.156089,'amu*angstrom^2'), symmetry=1, barrier=(3.5888,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.156089,'amu*angstrom^2'), symmetry=1, barrier=(3.5888,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.156089,'amu*angstrom^2'), symmetry=1, barrier=(3.5888,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.156089,'amu*angstrom^2'), symmetry=1, barrier=(3.5888,'kJ/mol'), semiclassical=False),
    ],
    spinMultiplicity = 5,
    opticalIsomers = 1,
    molecularWeight = (94.1112,'amu'),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[0.900944,0.0670114,-6.96506e-05,3.59982e-08,-6.02818e-12,81916.1,27], Tmin=(100,'K'), Tmax=(834.885,'K')), NASAPolynomial(coeffs=[12.4907,0.0219008,-7.31744e-06,1.16859e-09,-7.36356e-14,79617.9,-28.9963], Tmin=(834.885,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(680.151,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(286.849,'J/(mol*K)'), comment="""Thermo group additivity estimation: group(O2s-CsH) + group(Cs-CsCsCsH) + group(Cs-CtCsHH) + group(Cs-CsOsHH) + group(Cs-CsHHH) + group(Ct-CtCs) + group(Ct-CtH) + radical(Sec_Propargyl) + radical(CCJ2_triplet) + radical(CCJ(C)CO)"""),
)

species(
    label = '[CH]C=CC([CH])=CO(29851)',
    structure = SMILES('[CH]C=CC([CH])=CO'),
    E0 = (552.881,'kJ/mol'),
    modes = [
        HarmonicOscillator(frequencies=([2995,3010,3025,975,987.5,1000,1300,1337.5,1375,400,450,500,1630,1655,1680,3615,1277.5,1000,350,440,435,1725,200,800,960,1120,1280,1440,1600],'cm^-1')),
        HinderedRotor(inertia=(0.156089,'amu*angstrom^2'), symmetry=1, barrier=(3.5888,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.156089,'amu*angstrom^2'), symmetry=1, barrier=(3.5888,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.156089,'amu*angstrom^2'), symmetry=1, barrier=(3.5888,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.156089,'amu*angstrom^2'), symmetry=1, barrier=(3.5888,'kJ/mol'), semiclassical=False),
    ],
    spinMultiplicity = 5,
    opticalIsomers = 1,
    molecularWeight = (94.1112,'amu'),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[0.169327,0.0679276,-1.67006e-05,-4.06762e-08,2.48813e-11,66649.3,26.2291], Tmin=(100,'K'), Tmax=(919.275,'K')), NASAPolynomial(coeffs=[21.5074,0.0176523,-4.13097e-06,5.85179e-10,-4.00673e-14,60927.4,-84.7025], Tmin=(919.275,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(552.881,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(291.007,'J/(mol*K)'), comment="""Thermo group additivity estimation: group(O2s-(Cds-Cd)H) + group(Cs-(Cds-Cds)HHH) + group(Cs-(Cds-Cds)HHH) + group(Cds-Cds(Cds-Cds)Cs) + group(Cds-CdsCsH) + group(Cds-Cds(Cds-Cds)H) + group(Cds-CdsOsH) + radical(AllylJ2_triplet) + radical(AllylJ2_triplet)"""),
)

species(
    label = '[CH][C]1C=C([CH]O)C1(29852)',
    structure = SMILES('[CH][C]1C=C([CH]O)C1'),
    E0 = (618.016,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
    molecularWeight = (94.1112,'amu'),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[1.69557,0.030973,4.88649e-05,-9.24713e-08,3.95363e-11,74431.1,23.937], Tmin=(100,'K'), Tmax=(942.186,'K')), NASAPolynomial(coeffs=[17.4429,0.0139404,-3.33708e-06,5.89384e-10,-4.83744e-14,69252.4,-62.8323], Tmin=(942.186,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(618.016,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(295.164,'J/(mol*K)'), comment="""Thermo group additivity estimation: group(O2s-CsH) + group(Cs-(Cds-Cds)CsCsH) + group(Cs-(Cds-Cds)CsHH) + group(Cs-(Cds-Cds)OsHH) + group(Cs-CsHHH) + group(Cds-CdsCsCs) + group(Cds-CdsCsH) + ring(Cyclobutene) + radical(CCJ2_triplet) + radical(Allyl_T) + radical(C=CCJO)"""),
)

species(
    label = '[CH]=[C][CH]C(O)[C]=C(22023)',
    structure = SMILES('[CH]=[C][CH]C(O)[C]=C'),
    E0 = (683.816,'kJ/mol'),
    modes = [
        HarmonicOscillator(frequencies=([3120,650,792.5,1650,3615,1277.5,1000,3025,407.5,1350,352.5,2950,3100,1380,975,1025,1650,1380,1390,370,380,2900,435,1670,1700,300,440,244.921,709.027],'cm^-1')),
        HinderedRotor(inertia=(0.0589693,'amu*angstrom^2'), symmetry=1, barrier=(21.1138,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.0589364,'amu*angstrom^2'), symmetry=1, barrier=(21.1166,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.0592513,'amu*angstrom^2'), symmetry=1, barrier=(21.1179,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.919193,'amu*angstrom^2'), symmetry=1, barrier=(21.1341,'kJ/mol'), semiclassical=False),
    ],
    spinMultiplicity = 5,
    opticalIsomers = 1,
    molecularWeight = (94.1112,'amu'),
    collisionModel = TransportData(shapeIndex=2, epsilon=(3797.5,'J/mol'), sigma=(6.3848,'angstroms'), dipoleMoment=(0,'C*m'), polarizability=(0,'angstroms^3'), rotrelaxcollnum=0, comment="""Epsilon & sigma estimated with Tc=593.16 K, Pc=33.11 bar (from Joback method)"""),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[0.92828,0.057037,-3.30752e-05,-2.83369e-09,6.25036e-12,82364,31.4828], Tmin=(100,'K'), Tmax=(1009.59,'K')), NASAPolynomial(coeffs=[16.4734,0.0173712,-6.71531e-06,1.26956e-09,-9.20185e-14,78107.8,-49.1956], Tmin=(1009.59,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(683.816,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(291.007,'J/(mol*K)'), comment="""Thermo group additivity estimation: group(O2s-CsH) + group(Cs-(Cds-Cds)CsOsH) + group(Cs-(Cds-Cds)CsHH) + group(Cds-CdsCsH) + group(Cds-CdsCsH) + group(Cds-CdsHH) + group(Cds-CdsHH) + radical(Cds_P) + radical(C=CCJC(O)C=C) + radical(Cds_S) + radical(Cds_S)"""),
)

species(
    label = '[CH]C(=CO)C[C]=[CH](28754)',
    structure = SMILES('[CH]C(=CO)C[C]=[CH]'),
    E0 = (696.524,'kJ/mol'),
    modes = [
        HarmonicOscillator(frequencies=([3120,650,792.5,1650,3615,1277.5,1000,350,440,435,1725,2750,2850,1437.5,1250,1305,750,350,1685,370,3010,987.5,1337.5,450,1655,180,180,180,180],'cm^-1')),
        HinderedRotor(inertia=(2.04699,'amu*angstrom^2'), symmetry=1, barrier=(47.0643,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(2.04695,'amu*angstrom^2'), symmetry=1, barrier=(47.0634,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(2.04673,'amu*angstrom^2'), symmetry=1, barrier=(47.0583,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(2.04711,'amu*angstrom^2'), symmetry=1, barrier=(47.067,'kJ/mol'), semiclassical=False),
    ],
    spinMultiplicity = 5,
    opticalIsomers = 1,
    molecularWeight = (94.1112,'amu'),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[-0.399971,0.0793833,-7.78988e-05,3.77149e-08,-6.94236e-12,83945.8,31.2269], Tmin=(100,'K'), Tmax=(1488.65,'K')), NASAPolynomial(coeffs=[20.7528,0.0139206,-3.24576e-06,3.90526e-10,-2.05507e-14,78603.7,-76.0288], Tmin=(1488.65,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(696.524,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(291.007,'J/(mol*K)'), comment="""Thermo group additivity estimation: group(O2s-(Cds-Cd)H) + group(Cs-(Cds-Cds)(Cds-Cds)HH) + group(Cs-(Cds-Cds)HHH) + group(Cds-CdsCsCs) + group(Cds-CdsCsH) + group(Cds-CdsOsH) + group(Cds-CdsHH) + radical(Cds_S) + radical(AllylJ2_triplet) + radical(Cds_P)"""),
)

species(
    label = '[CH]C(C=[C][CH2])=CO(24296)',
    structure = SMILES('[CH]C(C=[C][CH2])=CO'),
    E0 = (538.094,'kJ/mol'),
    modes = [
        HarmonicOscillator(frequencies=([3000,3100,440,815,1455,1000,1685,370,2995,3025,975,1000,1300,1375,400,500,1630,1680,3615,1277.5,1000,350,440,435,1725,180,180,180,180],'cm^-1')),
        HinderedRotor(inertia=(2.07941,'amu*angstrom^2'), symmetry=1, barrier=(47.8096,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(2.07942,'amu*angstrom^2'), symmetry=1, barrier=(47.8099,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(2.0794,'amu*angstrom^2'), symmetry=1, barrier=(47.8095,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(2.07939,'amu*angstrom^2'), symmetry=1, barrier=(47.8092,'kJ/mol'), semiclassical=False),
    ],
    spinMultiplicity = 5,
    opticalIsomers = 1,
    molecularWeight = (94.1112,'amu'),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[0.416782,0.0636857,-1.85098e-05,-3.70321e-08,2.38442e-11,64861,26.2645], Tmin=(100,'K'), Tmax=(907.653,'K')), NASAPolynomial(coeffs=[21.1917,0.0128555,-1.80835e-06,1.33117e-10,-7.89505e-15,59412.2,-81.1903], Tmin=(907.653,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(538.094,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(291.007,'J/(mol*K)'), comment="""Thermo group additivity estimation: group(O2s-(Cds-Cd)H) + group(Cs-(Cds-Cds)HHH) + group(Cs-(Cds-Cds)HHH) + group(Cds-Cds(Cds-Cds)Cs) + group(Cds-CdsCsH) + group(Cds-Cds(Cds-Cds)H) + group(Cds-CdsOsH) + radical(C=CC=CCJ) + radical(Cds_S) + radical(AllylJ2_triplet)"""),
)

species(
    label = 'O[CH][C]1[CH][C]=CC1(29853)',
    structure = SMILES('O[CH][C]1[CH][C]=CC1'),
    E0 = (547.032,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
    molecularWeight = (94.1112,'amu'),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[1.75499,0.0368161,1.46351e-05,-4.43028e-08,1.91106e-11,65884.7,22.5027], Tmin=(100,'K'), Tmax=(991.543,'K')), NASAPolynomial(coeffs=[12.7563,0.0222742,-8.50593e-06,1.60608e-09,-1.16669e-13,62236.2,-37.8755], Tmin=(991.543,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(547.032,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(299.321,'J/(mol*K)'), comment="""Thermo group additivity estimation: group(O2s-CsH) + group(Cs-CsCsCsH) + group(Cs-(Cds-Cds)CsHH) + group(Cs-(Cds-Cds)CsHH) + group(Cs-CsOsHH) + group(Cds-CdsCsH) + group(Cds-CdsCsH) + ring(Cyclopentene) + radical(CCsJOH) + radical(cyclopentene-vinyl) + radical(CCJ(C)CO) + radical(cyclopentene-allyl)"""),
)

species(
    label = '[CH]=[C]CC(=C)C=O(19965)',
    structure = SMILES('[CH]=[C]CC(=C)C=O'),
    E0 = (442.882,'kJ/mol'),
    modes = [
        HarmonicOscillator(frequencies=([3120,650,792.5,1650,2950,3100,1380,975,1025,1650,350,440,435,1725,2750,2850,1437.5,1250,1305,750,350,2782.5,750,1395,475,1775,1000,1685,370,310.195],'cm^-1')),
        HinderedRotor(inertia=(0.182925,'amu*angstrom^2'), symmetry=1, barrier=(12.488,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.182881,'amu*angstrom^2'), symmetry=1, barrier=(12.488,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.182856,'amu*angstrom^2'), symmetry=1, barrier=(12.4879,'kJ/mol'), semiclassical=False),
    ],
    spinMultiplicity = 3,
    opticalIsomers = 1,
    molecularWeight = (94.1112,'amu'),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[1.68419,0.0569697,-5.07712e-05,2.53407e-08,-5.58423e-12,53344.7,24.6538], Tmin=(100,'K'), Tmax=(1027.47,'K')), NASAPolynomial(coeffs=[7.63504,0.0338027,-1.69497e-05,3.39579e-09,-2.44662e-13,52121.9,-4.2169], Tmin=(1027.47,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(442.882,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(295.164,'J/(mol*K)'), comment="""Thermo group additivity estimation: group(Cs-(Cds-Cds)(Cds-Cds)HH) + group(Cd-CdCs(CO)) + group(Cds-CdsCsH) + group(Cds-O2d(Cds-Cds)H) + group(Cds-CdsHH) + group(Cds-CdsHH) + radical(Cds_S) + radical(Cds_P)"""),
)

species(
    label = '[CH]=[C]C1C(=C)C1O(29824)',
    structure = SMILES('[CH]=[C]C1C(=C)C1O'),
    E0 = (570.777,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
    molecularWeight = (94.1112,'amu'),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[0.843179,0.0628792,-5.86119e-05,2.78836e-08,-5.24512e-12,68767.8,23.9063], Tmin=(100,'K'), Tmax=(1289.72,'K')), NASAPolynomial(coeffs=[15.117,0.0186093,-7.12397e-06,1.269e-09,-8.61043e-14,65086,-48.5884], Tmin=(1289.72,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(570.777,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(299.321,'J/(mol*K)'), comment="""Thermo group additivity estimation: group(O2s-CsH) + group(Cs-(Cds-Cds)(Cds-Cds)CsH) + group(Cs-(Cds-Cds)CsOsH) + group(Cds-CdsCsCs) + group(Cds-CdsCsH) + group(Cds-CdsHH) + group(Cds-CdsHH) + ring(Methylene_cyclopropane) + radical(Cds_S) + radical(Cds_P)"""),
)

species(
    label = '[C]=[C]CC([CH2])=CO(28755)',
    structure = SMILES('[C]=[C]CC([CH2])=CO'),
    E0 = (788.345,'kJ/mol'),
    modes = [
        HarmonicOscillator(frequencies=([3615,1277.5,1000,350,440,435,1725,2750,2850,1437.5,1250,1305,750,350,3000,3100,440,815,1455,1000,1685,370,3010,987.5,1337.5,450,1655,180,180],'cm^-1')),
        HinderedRotor(inertia=(0.888868,'amu*angstrom^2'), symmetry=1, barrier=(20.4368,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.88878,'amu*angstrom^2'), symmetry=1, barrier=(20.4348,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.889128,'amu*angstrom^2'), symmetry=1, barrier=(20.4428,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.888952,'amu*angstrom^2'), symmetry=1, barrier=(20.4388,'kJ/mol'), semiclassical=False),
    ],
    spinMultiplicity = 5,
    opticalIsomers = 1,
    molecularWeight = (94.1112,'amu'),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[-0.198936,0.0779757,-8.42149e-05,4.36259e-08,-8.5128e-12,94979.5,29.6455], Tmin=(100,'K'), Tmax=(1393.22,'K')), NASAPolynomial(coeffs=[21.3137,0.00863934,-1.41145e-06,1.02434e-10,-2.88199e-15,89720.1,-78.6369], Tmin=(1393.22,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(788.345,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(291.007,'J/(mol*K)'), comment="""Thermo group additivity estimation: group(O2s-(Cds-Cd)H) + group(Cs-(Cds-Cds)(Cds-Cds)HH) + group(Cs-(Cds-Cds)HHH) + group(Cds-CdsCsCs) + group(Cds-CdsCsH) + group(Cds-CdsOsH) + group(Cds-CdsHH) + radical(Cds_S) + radical(Allyl_P) + radical(CdCdJ2_triplet)"""),
)

species(
    label = '[C]#C[CH]C([CH2])[CH]O(29854)',
    structure = SMILES('[C]#C[CH]C([CH2])[CH]O'),
    E0 = (801.886,'kJ/mol'),
    modes = [
        HarmonicOscillator(frequencies=([3000,3100,440,815,1455,1000,3000,3050,390,425,1340,1360,335,370,3615,1277.5,1000,2175,525,1380,1390,370,380,2900,435,200,800,1600],'cm^-1')),
        HinderedRotor(inertia=(0.156089,'amu*angstrom^2'), symmetry=1, barrier=(3.5888,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.156089,'amu*angstrom^2'), symmetry=1, barrier=(3.5888,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.156089,'amu*angstrom^2'), symmetry=1, barrier=(3.5888,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.156089,'amu*angstrom^2'), symmetry=1, barrier=(3.5888,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.156089,'amu*angstrom^2'), symmetry=1, barrier=(3.5888,'kJ/mol'), semiclassical=False),
    ],
    spinMultiplicity = 5,
    opticalIsomers = 1,
    molecularWeight = (94.1112,'amu'),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[0.337042,0.0860898,-0.000133794,1.0667e-07,-3.17236e-11,96571.5,29.888], Tmin=(100,'K'), Tmax=(1002.81,'K')), NASAPolynomial(coeffs=[10.9363,0.0224286,-6.58517e-06,8.39453e-10,-3.92204e-14,95520.9,-15.9165], Tmin=(1002.81,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(801.886,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(286.849,'J/(mol*K)'), comment="""Thermo group additivity estimation: group(O2s-CsH) + group(Cs-CsCsCsH) + group(Cs-CtCsHH) + group(Cs-CsOsHH) + group(Cs-CsHHH) + group(Ct-CtCs) + group(Ct-CtH) + radical(Acetyl) + radical(CCsJOH) + radical(Sec_Propargyl) + radical(Isobutyl)"""),
)

species(
    label = '[C]=[C]C=C([CH2])CO(29855)',
    structure = SMILES('[C]=[C]C=C([CH2])CO'),
    E0 = (740.858,'kJ/mol'),
    modes = [
        HarmonicOscillator(frequencies=([3615,1277.5,1000,350,440,435,1725,2750,2850,1437.5,1250,1305,750,350,3000,3100,440,815,1455,1000,1685,370,3010,987.5,1337.5,450,1655,180,2729.58],'cm^-1')),
        HinderedRotor(inertia=(0.658025,'amu*angstrom^2'), symmetry=1, barrier=(15.1293,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(3.15387,'amu*angstrom^2'), symmetry=1, barrier=(72.5136,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(3.15478,'amu*angstrom^2'), symmetry=1, barrier=(72.5347,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(3.15462,'amu*angstrom^2'), symmetry=1, barrier=(72.5309,'kJ/mol'), semiclassical=False),
    ],
    spinMultiplicity = 5,
    opticalIsomers = 1,
    molecularWeight = (94.1112,'amu'),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[0.863644,0.0703874,-8.49058e-05,5.77095e-08,-1.5819e-11,89216.3,27.1055], Tmin=(100,'K'), Tmax=(890.539,'K')), NASAPolynomial(coeffs=[10.7405,0.0260252,-1.01856e-05,1.77484e-09,-1.17007e-13,87457.1,-19.3999], Tmin=(890.539,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(740.858,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(291.007,'J/(mol*K)'), comment="""Thermo group additivity estimation: group(O2s-CsH) + group(Cs-(Cds-Cds)OsHH) + group(Cs-(Cds-Cds)HHH) + group(Cds-CdsCsCs) + group(Cds-Cds(Cds-Cds)H) + group(Cds-Cds(Cds-Cds)H) + group(Cds-CdsHH) + radical(CdCdJ2_triplet) + radical(C=CC=CCJ) + radical(C=CJC=C)"""),
)

species(
    label = '[C][C]=CC(C)=CO(29856)',
    structure = SMILES('[C][C]=CC(C)=CO'),
    E0 = (718.832,'kJ/mol'),
    modes = [
        HarmonicOscillator(frequencies=([1685,370,2750,2800,2850,1350,1500,750,1050,1375,1000,2995,3025,975,1000,1300,1375,400,500,1630,1680,3615,1277.5,1000,350,440,435,1725,180,180],'cm^-1')),
        HinderedRotor(inertia=(1.02782,'amu*angstrom^2'), symmetry=1, barrier=(23.6316,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(1.02816,'amu*angstrom^2'), symmetry=1, barrier=(23.6395,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(1.03299,'amu*angstrom^2'), symmetry=1, barrier=(23.7505,'kJ/mol'), semiclassical=False),
    ],
    spinMultiplicity = 5,
    opticalIsomers = 1,
    molecularWeight = (94.1112,'amu'),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[-0.680531,0.0837937,-9.38743e-05,4.9089e-08,-9.4771e-12,86640.5,25.9891], Tmin=(100,'K'), Tmax=(1466.41,'K')), NASAPolynomial(coeffs=[23.5677,0.00419027,1.32111e-06,-4.4877e-10,3.54739e-14,80976.1,-95.3429], Tmin=(1466.41,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(718.832,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(295.164,'J/(mol*K)'), comment="""Thermo group additivity estimation: group(O2s-(Cds-Cd)H) + group(Cs-(Cds-Cds)HHH) + group(Cs-(Cds-Cds)HHH) + group(Cds-Cds(Cds-Cds)Cs) + group(Cds-CdsCsH) + group(Cds-Cds(Cds-Cds)H) + group(Cds-CdsOsH) + radical(Cds_S) + radical(CJ3)"""),
)

species(
    label = '[C]#C(5143)',
    structure = SMILES('[C]#C'),
    E0 = (551.936,'kJ/mol'),
    modes = [
        HarmonicOscillator(frequencies=([750,770,3400,2100],'cm^-1')),
    ],
    spinMultiplicity = 2,
    opticalIsomers = 1,
    molecularWeight = (25.0293,'amu'),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[3.03853,0.0115449,-2.13263e-05,1.81932e-08,-5.41587e-12,66398,5.96675], Tmin=(100,'K'), Tmax=(1076.58,'K')), NASAPolynomial(coeffs=[4.00845,0.00206817,6.04935e-08,-1.17707e-10,1.2928e-14,66529.5,2.79656], Tmin=(1076.58,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(551.936,'kJ/mol'), Cp0=(29.1007,'J/(mol*K)'), CpInf=(62.3585,'J/(mol*K)'), comment="""Thermo group additivity estimation: group(Ct-CtH) + group(Ct-CtH) + radical(Acetyl)"""),
)

species(
    label = '[CH][C]([CH2])[CH]O(29857)',
    structure = SMILES('[CH][C]([CH2])[CH]O'),
    E0 = (672.996,'kJ/mol'),
    modes = [
        HarmonicOscillator(frequencies=([3000,3100,440,815,1455,1000,360,370,350,3615,1277.5,1000,3025,407.5,1350,352.5,494.691,1862.29,1862.31,1862.33],'cm^-1')),
        HinderedRotor(inertia=(0.398767,'amu*angstrom^2'), symmetry=1, barrier=(69.2217,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.0639438,'amu*angstrom^2'), symmetry=1, barrier=(11.111,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.0639924,'amu*angstrom^2'), symmetry=1, barrier=(11.1134,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.0281251,'amu*angstrom^2'), symmetry=1, barrier=(69.218,'kJ/mol'), semiclassical=False),
    ],
    spinMultiplicity = 6,
    opticalIsomers = 1,
    molecularWeight = (69.0819,'amu'),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[1.7936,0.0496403,-5.63457e-05,3.25203e-08,-6.33553e-12,81021.3,22.487], Tmin=(100,'K'), Tmax=(763.784,'K')), NASAPolynomial(coeffs=[9.42825,0.0168591,-6.11067e-06,1.01857e-09,-6.54351e-14,79645,-13.6636], Tmin=(763.784,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(672.996,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(216.176,'J/(mol*K)'), comment="""Thermo group additivity estimation: group(O2s-CsH) + group(Cs-CsCsCsH) + group(Cs-CsOsHH) + group(Cs-CsHHH) + group(Cs-CsHHH) + radical(Isobutyl) + radical(CCJ2_triplet) + radical(CCsJOH) + radical(CCJ(C)CO)"""),
)

species(
    label = '[CH]=C=[C]C([CH2])[CH]O(29858)',
    structure = SMILES('[CH]=C=[C]C([CH2])[CH]O'),
    E0 = (710.176,'kJ/mol'),
    modes = [
        HarmonicOscillator(frequencies=([3120,650,792.5,1650,3615,1277.5,1000,3025,407.5,1350,352.5,540,610,2055,3000,3100,440,815,1455,1000,1685,370,1380,1390,370,380,2900,435,282.882],'cm^-1')),
        HinderedRotor(inertia=(0.132863,'amu*angstrom^2'), symmetry=1, barrier=(7.54467,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.132862,'amu*angstrom^2'), symmetry=1, barrier=(7.54466,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.132863,'amu*angstrom^2'), symmetry=1, barrier=(7.54467,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(1.1594,'amu*angstrom^2'), symmetry=1, barrier=(65.8363,'kJ/mol'), semiclassical=False),
    ],
    spinMultiplicity = 5,
    opticalIsomers = 1,
    molecularWeight = (94.1112,'amu'),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[0.385328,0.083582,-0.000124095,9.65408e-08,-2.87883e-11,85540.8,31.0447], Tmin=(100,'K'), Tmax=(950.926,'K')), NASAPolynomial(coeffs=[11.8618,0.0226877,-8.13373e-06,1.28815e-09,-7.7311e-14,83928.7,-20.745], Tmin=(950.926,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(710.176,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(291.007,'J/(mol*K)'), comment="""Thermo group additivity estimation: group(O2s-CsH) + group(Cs-(Cds-Cdd-Cd)CsCsH) + group(Cs-CsOsHH) + group(Cs-CsHHH) + group(Cds-CdsCsH) + group(Cds-CdsHH) + group(Cdd-CdsCds) + radical(CCsJOH) + radical(Isobutyl) + radical(Cds_S) + radical(C=C=CJ)"""),
)

species(
    label = '[CH2]C#CC([CH2])=CO(28959)',
    structure = SMILES('[CH2]C#CC([CH2])=CO'),
    E0 = (287.184,'kJ/mol'),
    modes = [
        HarmonicOscillator(frequencies=([3000,3033.33,3066.67,3100,415,465,780,850,1435,1475,900,1100,3010,987.5,1337.5,450,1655,3615,1277.5,1000,2100,2250,500,550,350,440,435,1725,219.105],'cm^-1')),
        HinderedRotor(inertia=(0.795362,'amu*angstrom^2'), symmetry=1, barrier=(27.382,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.797018,'amu*angstrom^2'), symmetry=1, barrier=(27.4175,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.800078,'amu*angstrom^2'), symmetry=1, barrier=(27.3794,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.800028,'amu*angstrom^2'), symmetry=1, barrier=(27.3622,'kJ/mol'), semiclassical=False),
    ],
    spinMultiplicity = 3,
    opticalIsomers = 1,
    molecularWeight = (94.1112,'amu'),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[0.759286,0.056069,-1.39305e-05,-3.71547e-08,2.28159e-11,34671,24.712], Tmin=(100,'K'), Tmax=(920.536,'K')), NASAPolynomial(coeffs=[20.9562,0.00804506,-4.28103e-07,-3.90985e-11,4.41362e-16,29269,-80.1992], Tmin=(920.536,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(287.184,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(291.007,'J/(mol*K)'), comment="""Thermo group additivity estimation: group(O2s-(Cds-Cd)H) + group(Cs-(Cds-Cds)HHH) + group(Cs-CtHHH) + group(Cds-CdsCtCs) + group(Cds-CdsOsH) + group(Ct-CtCs) + group(Ct-Ct(Cds-Cds)) + radical(Allyl_P) + radical(Propargyl)"""),
)

species(
    label = '[CH2][C][CH]O(28660)',
    structure = SMILES('[CH2][C][CH]O'),
    E0 = (553.408,'kJ/mol'),
    modes = [
        HarmonicOscillator(frequencies=([3025,407.5,1350,352.5,3615,1277.5,1000,3000,3100,440,815,1455,1000,180,2262.1],'cm^-1')),
        HinderedRotor(inertia=(0.0475559,'amu*angstrom^2'), symmetry=1, barrier=(11.7885,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.515056,'amu*angstrom^2'), symmetry=1, barrier=(11.8422,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.400962,'amu*angstrom^2'), symmetry=1, barrier=(46.0504,'kJ/mol'), semiclassical=False),
    ],
    spinMultiplicity = 5,
    opticalIsomers = 1,
    molecularWeight = (56.0633,'amu'),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[2.12807,0.0437773,-6.51486e-05,5.20447e-08,-1.62809e-11,66624.5,18.0235], Tmin=(100,'K'), Tmax=(877.804,'K')), NASAPolynomial(coeffs=[7.75761,0.0135713,-5.75184e-06,1.02548e-09,-6.75894e-14,65811.6,-7.40289], Tmin=(877.804,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(553.408,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(170.447,'J/(mol*K)'), comment="""Thermo group additivity estimation: group(O2s-CsH) + group(Cs-CsCsHH) + group(Cs-CsOsHH) + group(Cs-CsHHH) + radical(CCJ2_triplet) + radical(CCsJOH) + radical(RCCJ)"""),
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
    E0 = (547.147,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS2',
    E0 = (792.861,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS3',
    E0 = (668.28,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS4',
    E0 = (737.94,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS5',
    E0 = (772.454,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS6',
    E0 = (808.533,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS7',
    E0 = (767.483,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS8',
    E0 = (761.611,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS9',
    E0 = (776.125,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS10',
    E0 = (774.99,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS11',
    E0 = (769.473,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS12',
    E0 = (659.273,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS13',
    E0 = (580.188,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS14',
    E0 = (698.668,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS15',
    E0 = (924.805,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS16',
    E0 = (957.947,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS17',
    E0 = (751.624,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS18',
    E0 = (739.483,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS19',
    E0 = (736.8,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS20',
    E0 = (722.507,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS21',
    E0 = (610.547,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS22',
    E0 = (610.547,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS23',
    E0 = (572.121,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS24',
    E0 = (572.121,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS25',
    E0 = (554.679,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS26',
    E0 = (555.055,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS27',
    E0 = (555.432,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS28',
    E0 = (987.234,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS29',
    E0 = (1100.11,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS30',
    E0 = (1072.31,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS31',
    E0 = (1162.12,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS32',
    E0 = (1023.09,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS33',
    E0 = (1002.53,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS34',
    E0 = (1069.96,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS35',
    E0 = (867.709,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS36',
    E0 = (701.038,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS37',
    E0 = (802.371,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS38',
    E0 = (976.873,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS39',
    E0 = (1048.42,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS40',
    E0 = (666.927,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS41',
    E0 = (685.203,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS42',
    E0 = (975.613,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS43',
    E0 = (1010.48,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS44',
    E0 = (860.337,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS45',
    E0 = (691.014,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS46',
    E0 = (867.429,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS47',
    E0 = (732.523,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS48',
    E0 = (740.548,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS49',
    E0 = (579.501,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS50',
    E0 = (552.542,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS51',
    E0 = (825.336,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS52',
    E0 = (956.314,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS53',
    E0 = (685.203,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS54',
    E0 = (778.291,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS55',
    E0 = (888.653,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS56',
    E0 = (738.356,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS57',
    E0 = (579.501,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS58',
    E0 = (572.121,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS59',
    E0 = (570.879,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS60',
    E0 = (972.518,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS61',
    E0 = (954.36,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS62',
    E0 = (805.579,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS63',
    E0 = (760.404,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS64',
    E0 = (1259.25,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS65',
    E0 = (887.862,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS66',
    E0 = (610.547,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS67',
    E0 = (1079.41,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

reaction(
    label = 'reaction1',
    reactants = ['[CH]=[C]C=C([CH2])[CH]O(21903)'],
    products = ['C=C=CO(12571)', '[CH]=C=[CH](18734)'],
    transitionState = 'TS1',
    kinetics = Arrhenius(A=(5e+12,'s^-1'), n=0, Ea=(0,'kJ/mol'), T0=(1,'K'), Tmin=(300,'K'), Tmax=(1500,'K'), comment="""Exact match found for rate rule [RJJ]
Euclidian distance = 0
family: 1,4_Linear_birad_scission"""),
)

reaction(
    label = 'reaction2',
    reactants = ['[CH]=[C]C=C([CH2])[CH]O(21903)'],
    products = ['[CH2][C]([CH]O)C1[C]=C1(29830)'],
    transitionState = 'TS2',
    kinetics = Arrhenius(A=(2.54e+10,'s^-1'), n=0.69, Ea=(245.714,'kJ/mol'), T0=(1,'K'), comment="""Estimated using an average for rate rule [R4;doublebond_intra;radadd_intra_cdsingleH]
Euclidian distance = 0
family: Intra_R_Add_Exocyclic"""),
)

reaction(
    label = 'reaction3',
    reactants = ['H(8)', '[CH]=C=CC([CH2])=C[O](21884)'],
    products = ['[CH]=[C]C=C([CH2])[CH]O(21903)'],
    transitionState = 'TS3',
    kinetics = Arrhenius(A=(373000,'cm^3/(mol*s)'), n=2.53, Ea=(20.92,'kJ/mol'), T0=(1,'K'), Tmin=(300,'K'), Tmax=(1500,'K'), comment="""From training reaction 2834 used for Od_CO-CdH;HJ
Exact match found for rate rule [Od_CO-CdH;HJ]
Euclidian distance = 0
family: R_Addition_MultipleBond"""),
)

reaction(
    label = 'reaction4',
    reactants = ['H(8)', '[CH]=C=C=C([CH2])[CH]O(29831)'],
    products = ['[CH]=[C]C=C([CH2])[CH]O(21903)'],
    transitionState = 'TS4',
    kinetics = Arrhenius(A=(8.22e+08,'cm^3/(mol*s)'), n=1.533, Ea=(7.77387,'kJ/mol'), T0=(1,'K'), comment="""From training reaction 192 used for Cd_R;HJ
Exact match found for rate rule [Cd_R;HJ]
Euclidian distance = 0
family: R_Addition_MultipleBond"""),
)

reaction(
    label = 'reaction5',
    reactants = ['CH2(T)(28)', '[CH]=C=C[C]=CO(23383)'],
    products = ['[CH]=[C]C=C([CH2])[CH]O(21903)'],
    transitionState = 'TS5',
    kinetics = Arrhenius(A=(10.4286,'m^3/(mol*s)'), n=1.99344, Ea=(11.9048,'kJ/mol'), T0=(1,'K'), comment="""Estimated using template [Ca_Cds;YJ] for rate rule [Ca_Cds;CH2_triplet]
Euclidian distance = 2.0
family: R_Addition_MultipleBond"""),
)

reaction(
    label = 'reaction6',
    reactants = ['[CH]O(5471)', '[CH]=C=C[C]=C(19277)'],
    products = ['[CH]=[C]C=C([CH2])[CH]O(21903)'],
    transitionState = 'TS6',
    kinetics = Arrhenius(A=(1.70446,'m^3/(mol*s)'), n=2.07639, Ea=(14.6531,'kJ/mol'), T0=(1,'K'), comment="""Estimated using template [Ca_Cds-HH;YJ] for rate rule [Ca_Cds-HH;Y_1centerbirad]
Euclidian distance = 1.0
family: R_Addition_MultipleBond"""),
)

reaction(
    label = 'reaction7',
    reactants = ['[CH]=[C]C=C([CH2])C[O](22137)'],
    products = ['[CH]=[C]C=C([CH2])[CH]O(21903)'],
    transitionState = 'TS7',
    kinetics = Arrhenius(A=(2.42705e+08,'s^-1'), n=1.50595, Ea=(111.926,'kJ/mol'), T0=(1,'K'), comment="""Estimated using template [R2H_S;Y_rad_out;Cs_H_out_H/Cd] for rate rule [R2H_S;O_rad_out;Cs_H_out_H/Cd]
Euclidian distance = 1.0
Multiplied by reaction path degeneracy 2.0
family: intra_H_migration"""),
)

reaction(
    label = 'reaction8',
    reactants = ['[CH]=[C]C=C([CH2])[CH]O(21903)'],
    products = ['[CH]C=[C]C([CH2])=CO(29832)'],
    transitionState = 'TS8',
    kinetics = Arrhenius(A=(1.806e+09,'s^-1'), n=1.172, Ea=(214.463,'kJ/mol'), T0=(1,'K'), comment="""Estimated using an average for rate rule [R2H_S;Cd_rad_out;Cd_H_out_doubleC]
Euclidian distance = 0
family: intra_H_migration"""),
)

reaction(
    label = 'reaction9',
    reactants = ['[CH]=[C][C]=C([CH2])CO(29833)'],
    products = ['[CH]=[C]C=C([CH2])[CH]O(21903)'],
    transitionState = 'TS9',
    kinetics = Arrhenius(A=(1.47715e+10,'s^-1'), n=0.8, Ea=(147.277,'kJ/mol'), T0=(1,'K'), comment="""Estimated using template [R3H_DS;Cd_rad_out;Cs_H_out_1H] for rate rule [R3H_DS;Cd_rad_out;Cs_H_out_H/NonDeO]
Euclidian distance = 1.0
Multiplied by reaction path degeneracy 2.0
family: intra_H_migration"""),
)

reaction(
    label = 'reaction10',
    reactants = ['[CH][C]=[C]C(C)=CO(29834)'],
    products = ['[CH]=[C]C=C([CH2])[CH]O(21903)'],
    transitionState = 'TS10',
    kinetics = Arrhenius(A=(5.17353e+06,'s^-1'), n=1.89718, Ea=(155.956,'kJ/mol'), T0=(1,'K'), comment="""Estimated using an average for rate rule [R3H_DS;Cd_rad_out;Cs_H_out_2H]
Euclidian distance = 0
Multiplied by reaction path degeneracy 3.0
family: intra_H_migration"""),
)

reaction(
    label = 'reaction11',
    reactants = ['[CH]=[C]C=C([CH2])[CH]O(21903)'],
    products = ['[CH2]C(=[C][C]=C)[CH]O(28964)'],
    transitionState = 'TS11',
    kinetics = Arrhenius(A=(191.5,'s^-1'), n=3.05, Ea=(222.325,'kJ/mol'), T0=(1,'K'), comment="""Estimated using template [R3Hall;Cd_rad_out_singleH;Cd_H_out_doubleC] for rate rule [R3HJ;Cd_rad_out_singleH;Cd_H_out_doubleC]
Euclidian distance = 1.0
family: intra_H_migration"""),
)

reaction(
    label = 'reaction12',
    reactants = ['[CH][C]=C[C](C)C=O(21886)'],
    products = ['[CH]=[C]C=C([CH2])[CH]O(21903)'],
    transitionState = 'TS12',
    kinetics = Arrhenius(A=(2.3e+09,'s^-1','*|/',2.51189), n=0.75, Ea=(97.0688,'kJ/mol'), T0=(1,'K'), comment="""Estimated using template [R4Hall;O_rad_out;Cs_H_out_2H] for rate rule [R4HJ_1;O_rad_out;Cs_H_out_2H]
Euclidian distance = 1.0
Multiplied by reaction path degeneracy 3.0
family: intra_H_migration"""),
)

reaction(
    label = 'reaction13',
    reactants = ['[CH]=[C]C=C([CH2])[CH]O(21903)'],
    products = ['[CH]C=CC([CH2])=C[O](21890)'],
    transitionState = 'TS13',
    kinetics = Arrhenius(A=(136000,'s^-1'), n=1.9199, Ea=(33.0402,'kJ/mol'), T0=(1,'K'), comment="""Estimated using template [R5Hall;Cd_rad_out;XH_out] for rate rule [R5HJ_3;Cd_rad_out;O_H_out]
Euclidian distance = 1.41421356237
family: intra_H_migration"""),
)

reaction(
    label = 'reaction14',
    reactants = ['[CH]=[C]C=C([CH2])[CH]O(21903)'],
    products = ['[CH2][C]=CC([CH2])=C[O](20121)'],
    transitionState = 'TS14',
    kinetics = Arrhenius(A=(1.86943e+06,'s^-1'), n=1.85754, Ea=(151.521,'kJ/mol'), T0=(1,'K'), comment="""Estimated using template [RnH;Cd_rad_out_singleH;XH_out] for rate rule [R6Hall;Cd_rad_out_singleH;O_H_out]
Euclidian distance = 1.41421356237
family: intra_H_migration"""),
)

reaction(
    label = 'reaction15',
    reactants = ['H(8)', '[CH][C]=CC([CH2])=C[O](21892)'],
    products = ['[CH]=[C]C=C([CH2])[CH]O(21903)'],
    transitionState = 'TS15',
    kinetics = Arrhenius(A=(4.34078e+06,'m^3/(mol*s)'), n=0.278577, Ea=(0,'kJ/mol'), T0=(1,'K'), comment="""Estimated using an average for rate rule [Y_rad;H_rad]
Euclidian distance = 0
family: R_Recombination
Ea raised from -1.4 to 0 kJ/mol."""),
)

reaction(
    label = 'reaction16',
    reactants = ['H(8)', '[CH]=[C][C]=C([CH2])[CH]O(29835)'],
    products = ['[CH]=[C]C=C([CH2])[CH]O(21903)'],
    transitionState = 'TS16',
    kinetics = Arrhenius(A=(4.34078e+06,'m^3/(mol*s)'), n=0.278577, Ea=(0,'kJ/mol'), T0=(1,'K'), comment="""Estimated using an average for rate rule [Y_rad;H_rad]
Euclidian distance = 0
family: R_Recombination
Ea raised from -1.4 to 0 kJ/mol."""),
)

reaction(
    label = 'reaction17',
    reactants = ['[CH]=[C]C=C([CH2])[CH]O(21903)'],
    products = ['[CH]=[C]C1[C]([CH2])C1O(29836)'],
    transitionState = 'TS17',
    kinetics = Arrhenius(A=(4.00692e+11,'s^-1'), n=0.347401, Ea=(204.477,'kJ/mol'), T0=(1,'K'), comment="""Estimated using template [R3_D;doublebond_intra;radadd_intra_cs] for rate rule [R3_D;doublebond_intra;radadd_intra_csHO]
Euclidian distance = 2.0
family: Intra_R_Add_Endocyclic
Ea raised from 202.4 to 204.5 kJ/mol to match endothermicity of reaction."""),
)

reaction(
    label = 'reaction18',
    reactants = ['[CH]=[C]C=C([CH2])[CH]O(21903)'],
    products = ['[CH]=[C]C1C[C]1[CH]O(29837)'],
    transitionState = 'TS18',
    kinetics = Arrhenius(A=(2.03537e+13,'s^-1'), n=-0.296394, Ea=(192.336,'kJ/mol'), T0=(1,'K'), Tmin=(303.03,'K'), Tmax=(2000,'K'), comment="""From training reaction 135 used for R3_D;doublebond_intra;radadd_intra_cs2H
Exact match found for rate rule [R3_D;doublebond_intra;radadd_intra_cs2H]
Euclidian distance = 0
family: Intra_R_Add_Endocyclic
Ea raised from 191.4 to 192.3 kJ/mol to match endothermicity of reaction."""),
)

reaction(
    label = 'reaction19',
    reactants = ['[CH]=[C]C=C([CH2])[CH]O(21903)'],
    products = ['[CH]=C1[CH]C1([CH2])[CH]O(29797)'],
    transitionState = 'TS19',
    kinetics = Arrhenius(A=(2.42595e+10,'s^-1'), n=0.7335, Ea=(189.653,'kJ/mol'), T0=(1,'K'), comment="""Estimated using an average for rate rule [R3_D;doublebond_intra_pri;radadd_intra]
Euclidian distance = 0
family: Intra_R_Add_Endocyclic"""),
)

reaction(
    label = 'reaction20',
    reactants = ['[CH]=[C]C=C([CH2])[CH]O(21903)'],
    products = ['[CH2]C1([CH]O)[CH][C]=C1(29838)'],
    transitionState = 'TS20',
    kinetics = Arrhenius(A=(1.953e+11,'s^-1'), n=0.387, Ea=(175.36,'kJ/mol'), T0=(1,'K'), comment="""Estimated using an average for rate rule [R4_linear;doublebond_intra_pri;radadd_intra_cdsingleH]
Euclidian distance = 0
family: Intra_R_Add_Endocyclic
Ea raised from 174.6 to 175.4 kJ/mol to match endothermicity of reaction."""),
)

reaction(
    label = 'reaction21',
    reactants = ['[CH]=[C]C=C([CH2])[CH]O(21903)'],
    products = ['[CH]=C=C=C([CH2])CO(29839)'],
    transitionState = 'TS21',
    kinetics = Arrhenius(A=(7.437e+08,'s^-1'), n=1.045, Ea=(63.4002,'kJ/mol'), T0=(1,'K'), comment="""Estimated using an average for rate rule [R3radExo;Y_rad_NDe;XH_Rrad]
Euclidian distance = 0
family: Intra_Disproportionation"""),
)

reaction(
    label = 'reaction22',
    reactants = ['[CH]=[C]C=C([CH2])[CH]O(21903)'],
    products = ['[CH]=C=C=C(C)[CH]O(29840)'],
    transitionState = 'TS22',
    kinetics = Arrhenius(A=(7.437e+08,'s^-1'), n=1.045, Ea=(63.4002,'kJ/mol'), T0=(1,'K'), comment="""Estimated using an average for rate rule [R3radExo;Y_rad;XH_Rrad]
Euclidian distance = 0
family: Intra_Disproportionation"""),
)

reaction(
    label = 'reaction23',
    reactants = ['[CH]=[C]C=C([CH2])[CH]O(21903)'],
    products = ['[CH]=C=C[C](C)C=O(21909)'],
    transitionState = 'TS23',
    kinetics = Arrhenius(A=(5.14222e+08,'s^-1'), n=0.311, Ea=(24.9733,'kJ/mol'), T0=(1,'K'), comment="""Estimated using template [R4;Y_rad;XH_Rrad] for rate rule [R4radEndo;Y_rad;XH_Rrad]
Euclidian distance = 1.0
family: Intra_Disproportionation"""),
)

reaction(
    label = 'reaction24',
    reactants = ['[CH]=[C]C=C([CH2])[CH]O(21903)'],
    products = ['[CH]C=CC(=C)C=O(21897)'],
    transitionState = 'TS24',
    kinetics = Arrhenius(A=(2.1261e+09,'s^-1'), n=0.137, Ea=(24.9733,'kJ/mol'), T0=(1,'K'), comment="""Estimated using template [R5;Y_rad;XH_Rrad] for rate rule [R5radEndo;Y_rad;XH_Rrad]
Euclidian distance = 1.0
family: Intra_Disproportionation"""),
)

reaction(
    label = 'reaction25',
    reactants = ['[CH]=[C]C=C([CH2])[CH]O(21903)'],
    products = ['C#C[CH][C]1CC1O(29841)'],
    transitionState = 'TS25',
    kinetics = Arrhenius(A=(2e+12,'s^-1'), n=0, Ea=(7.5312,'kJ/mol'), T0=(1,'K'), comment="""Estimated using template [Rn;C_rad_out_1H;Cpri_rad_out_2H] for rate rule [R3_SS;C_rad_out_H/NonDeO;Cpri_rad_out_2H]
Euclidian distance = 2.2360679775
family: Birad_recombination"""),
)

reaction(
    label = 'reaction26',
    reactants = ['[CH]=[C]C=C([CH2])[CH]O(21903)'],
    products = ['[CH]=C1[CH]C(=C)C1O(29795)'],
    transitionState = 'TS26',
    kinetics = Arrhenius(A=(1.8e+12,'s^-1'), n=-0.1525, Ea=(7.90776,'kJ/mol'), T0=(1,'K'), comment="""Estimated using average of templates [Rn;C_rad_out_1H;Ypri_rad_out] + [R4;C_rad_out_single;Ypri_rad_out] for rate rule [R4_SDS;C_rad_out_H/NonDeO;Ypri_rad_out]
Euclidian distance = 2.2360679775
family: Birad_recombination"""),
)

reaction(
    label = 'reaction27',
    reactants = ['[CH]=[C]C=C([CH2])[CH]O(21903)'],
    products = ['[CH]=C1[CH]C(=CO)C1(29790)'],
    transitionState = 'TS27',
    kinetics = Arrhenius(A=(1.62e+12,'s^-1'), n=-0.305, Ea=(8.28432,'kJ/mol'), T0=(1,'K'), comment="""Estimated using template [R4;C_rad_out_2H;Ypri_rad_out] for rate rule [R4_SDS;C_rad_out_2H;Ypri_rad_out]
Euclidian distance = 1.0
family: Birad_recombination"""),
)

reaction(
    label = 'reaction28',
    reactants = ['OH(D)(132)', '[CH]C([CH2])=C[C]=[CH](19887)'],
    products = ['[CH]=[C]C=C([CH2])[CH]O(21903)'],
    transitionState = 'TS28',
    kinetics = Arrhenius(A=(43.5839,'m^3/(mol*s)'), n=1.88017, Ea=(5.1666,'kJ/mol'), T0=(1,'K'), comment="""Estimated using template [O_rad;Birad] for rate rule [O_pri_rad;Birad]
Euclidian distance = 1.0
family: Birad_R_Recombination"""),
)

reaction(
    label = 'reaction29',
    reactants = ['[CH]O(5471)', '[CH]=[C]C=[C][CH2](19274)'],
    products = ['[CH]=[C]C=C([CH2])[CH]O(21903)'],
    transitionState = 'TS29',
    kinetics = Arrhenius(A=(1.14854e+06,'m^3/(mol*s)'), n=0.575199, Ea=(34.3157,'kJ/mol'), T0=(1,'K'), comment="""Estimated using an average for rate rule [Y_rad;Birad]
Euclidian distance = 0
family: Birad_R_Recombination"""),
)

reaction(
    label = 'reaction30',
    reactants = ['CH2(T)(28)', '[CH][C]=C[C]=CO(29842)'],
    products = ['[CH]=[C]C=C([CH2])[CH]O(21903)'],
    transitionState = 'TS30',
    kinetics = Arrhenius(A=(1.14854e+06,'m^3/(mol*s)'), n=0.575199, Ea=(34.3157,'kJ/mol'), T0=(1,'K'), comment="""Estimated using an average for rate rule [Y_rad;Birad]
Euclidian distance = 0
family: Birad_R_Recombination"""),
)

reaction(
    label = 'reaction31',
    reactants = ['[CH]C([CH2])=CO(19098)', '[C]=[CH](18830)'],
    products = ['[CH]=[C]C=C([CH2])[CH]O(21903)'],
    transitionState = 'TS31',
    kinetics = Arrhenius(A=(1.14854e+06,'m^3/(mol*s)'), n=0.575199, Ea=(34.3157,'kJ/mol'), T0=(1,'K'), comment="""Estimated using template [Y_rad;Birad] for rate rule [Cd_pri_rad;Birad]
Euclidian distance = 2.0
family: Birad_R_Recombination"""),
)

reaction(
    label = 'reaction32',
    reactants = ['H(8)', '[CH][C]=CC([CH2])=[C]O(29843)'],
    products = ['[CH]=[C]C=C([CH2])[CH]O(21903)'],
    transitionState = 'TS32',
    kinetics = Arrhenius(A=(1e+07,'m^3/(mol*s)'), n=0, Ea=(0,'kJ/mol'), T0=(1,'K'), comment="""Estimated using an average for rate rule [H_rad;Birad]
Euclidian distance = 0
family: Birad_R_Recombination"""),
)

reaction(
    label = 'reaction33',
    reactants = ['H(8)', '[CH][C]=CC([CH])=CO(24780)'],
    products = ['[CH]=[C]C=C([CH2])[CH]O(21903)'],
    transitionState = 'TS33',
    kinetics = Arrhenius(A=(1e+07,'m^3/(mol*s)'), n=0, Ea=(0,'kJ/mol'), T0=(1,'K'), comment="""Estimated using an average for rate rule [H_rad;Birad]
Euclidian distance = 0
family: Birad_R_Recombination"""),
)

reaction(
    label = 'reaction34',
    reactants = ['H(8)', '[C]=[C]C=C([CH2])[CH]O(29844)'],
    products = ['[CH]=[C]C=C([CH2])[CH]O(21903)'],
    transitionState = 'TS34',
    kinetics = Arrhenius(A=(1e+07,'m^3/(mol*s)'), n=0, Ea=(0,'kJ/mol'), T0=(1,'K'), comment="""Estimated using an average for rate rule [H_rad;Birad]
Euclidian distance = 0
family: Birad_R_Recombination"""),
)

reaction(
    label = 'reaction35',
    reactants = ['C=C=CO(12571)', '[CH][C]=[CH](21256)'],
    products = ['[CH]=[C]C=C([CH2])[CH]O(21903)'],
    transitionState = 'TS35',
    kinetics = Arrhenius(A=(0.00173894,'m^3/(mol*s)'), n=2.67356, Ea=(32.0272,'kJ/mol'), T0=(1,'K'), comment="""Estimated using an average for rate rule [Ca_Cds-HH;CJ]
Euclidian distance = 0
Multiplied by reaction path degeneracy 2.0
family: R_Addition_MultipleBond"""),
)

reaction(
    label = 'reaction36',
    reactants = ['[CH2][C]=CO(18753)', '[CH]=C=[CH](18734)'],
    products = ['[CH]=[C]C=C([CH2])[CH]O(21903)'],
    transitionState = 'TS36',
    kinetics = Arrhenius(A=(1.04713,'m^3/(mol*s)'), n=2.10494, Ea=(22.6844,'kJ/mol'), T0=(1,'K'), comment="""Estimated using an average for rate rule [Ct_Ct;CJ]
Euclidian distance = 0
Multiplied by reaction path degeneracy 2.0
family: R_Addition_MultipleBond"""),
)

reaction(
    label = 'reaction37',
    reactants = ['[CH]=[C]C=C([CH2])[CH]O(21903)'],
    products = ['[CH][C]=CC(C)=[C]O(29845)'],
    transitionState = 'TS37',
    kinetics = Arrhenius(A=(1.91e+11,'s^-1'), n=0.63, Ea=(255.224,'kJ/mol'), T0=(1,'K'), Tmin=(300,'K'), Tmax=(1500,'K'), comment="""From training reaction 199 used for R3H_SD;C_rad_out_2H;Cd_H_out_singleNd
Exact match found for rate rule [R3H_SD;C_rad_out_2H;Cd_H_out_singleNd]
Euclidian distance = 0
family: intra_H_migration"""),
)

reaction(
    label = 'reaction38',
    reactants = ['[CH]C=CC([CH2])=[C]O(29846)'],
    products = ['[CH]=[C]C=C([CH2])[CH]O(21903)'],
    transitionState = 'TS38',
    kinetics = Arrhenius(A=(1.13764e+12,'s^-1'), n=1.09983, Ea=(403.433,'kJ/mol'), T0=(1,'K'), comment="""Estimated using template [R4H_DSD;Cd_rad_out_single;XH_out] for rate rule [R4H_DSD;Cd_rad_out_singleNd;XH_out]
Euclidian distance = 1.0
family: intra_H_migration"""),
)

reaction(
    label = 'reaction39',
    reactants = ['[CH2][C]=CO(18753)', '[CH][C]=[CH](21256)'],
    products = ['[CH]=[C]C=C([CH2])[CH]O(21903)'],
    transitionState = 'TS39',
    kinetics = Arrhenius(A=(3.9578e+07,'m^3/(mol*s)'), n=-0.126319, Ea=(0,'kJ/mol'), T0=(1,'K'), comment="""Estimated using an average for rate rule [Y_rad;Y_rad]
Euclidian distance = 0
Multiplied by reaction path degeneracy 2.0
family: R_Recombination
Ea raised from -15.6 to -15.6 kJ/mol.
Ea raised from -15.6 to 0 kJ/mol."""),
)

reaction(
    label = 'reaction40',
    reactants = ['[CH]=[C]C=C([CH2])[CH]O(21903)'],
    products = ['[CH][C]=C[C]1CC1O(29847)'],
    transitionState = 'TS40',
    kinetics = Arrhenius(A=(2.03537e+13,'s^-1'), n=-0.296394, Ea=(119.779,'kJ/mol'), T0=(1,'K'), Tmin=(303.03,'K'), Tmax=(2000,'K'), comment="""Estimated using template [R3_D;doublebond_intra;radadd_intra_cs2H] for rate rule [R3_D;doublebond_intra_secDe_HNd;radadd_intra_cs2H]
Euclidian distance = 2.0
family: Intra_R_Add_Endocyclic
Ea raised from 117.9 to 119.8 kJ/mol to match endothermicity of reaction."""),
)

reaction(
    label = 'reaction41',
    reactants = ['[CH]=[C]C=C([CH2])[CH]O(21903)'],
    products = ['[CH]C1=C[C]([CH2])C1O(29848)'],
    transitionState = 'TS41',
    kinetics = Arrhenius(A=(1.953e+11,'s^-1'), n=0.387, Ea=(138.055,'kJ/mol'), T0=(1,'K'), comment="""Estimated using an average for rate rule [R4_D_D;doublebond_intra;radadd_intra]
Euclidian distance = 0
family: Intra_R_Add_Endocyclic"""),
)

reaction(
    label = 'reaction42',
    reactants = ['[CH][C]=C[C]([CH2])C[O](21898)'],
    products = ['[CH]=[C]C=C([CH2])[CH]O(21903)'],
    transitionState = 'TS42',
    kinetics = Arrhenius(A=(2.94659e+10,'s^-1'), n=0.2847, Ea=(27.8529,'kJ/mol'), T0=(1,'K'), comment="""Estimated using average of templates [Rn;Y_rad;XH_Rrad_De] + [R2radExo;Y_rad;XH_Rrad] for rate rule [R2radExo;Y_rad;XH_Rrad_De]
Euclidian distance = 1.0
Multiplied by reaction path degeneracy 2.0
family: Intra_Disproportionation"""),
)

reaction(
    label = 'reaction43',
    reactants = ['[CH][C]=[C]C([CH2])[CH]O(21900)'],
    products = ['[CH]=[C]C=C([CH2])[CH]O(21903)'],
    transitionState = 'TS43',
    kinetics = Arrhenius(A=(1.949e+11,'s^-1'), n=0.486, Ea=(22.8614,'kJ/mol'), T0=(1,'K'), comment="""Estimated using an average for rate rule [R2radExo;Y_rad;XH_Rrad_NDe]
Euclidian distance = 0
family: Intra_Disproportionation"""),
)

reaction(
    label = 'reaction23',
    reactants = ['[CH][C]=CC[C]=CO(21983)'],
    products = ['[CH]=[C]C=C([CH2])[CH]O(21903)'],
    transitionState = 'TS44',
    kinetics = Arrhenius(A=(1.74842e+09,'s^-1'), n=1.084, Ea=(170.038,'kJ/mol'), T0=(1,'K'), comment="""Estimated using average of templates [cCsCJ;CdsJ;C] + [cCs(-HH)CJ;CJ;C] for rate rule [cCs(-HH)CJ;CdsJ;C]
Euclidian distance = 1.0
family: 1,2_shiftC"""),
)

reaction(
    label = 'reaction31',
    reactants = ['[CH]=[C]C=C([CH2])[CH]O(21903)'],
    products = ['[CH][C]=CC([CH2])C=O(19974)'],
    transitionState = 'TS45',
    kinetics = Arrhenius(A=(605.045,'s^-1'), n=2.96, Ea=(143.867,'kJ/mol'), T0=(1,'K'), comment="""Estimated using an average for rate rule [R_ROR;R1_doublebond;R2_doublebond_H;R_O_H]
Euclidian distance = 0
family: ketoenol"""),
)

reaction(
    label = 'reaction46',
    reactants = ['[CH]=[C]CC([CH2])=[C]O(28748)'],
    products = ['[CH]=[C]C=C([CH2])[CH]O(21903)'],
    transitionState = 'TS46',
    kinetics = Arrhenius(A=(9.9395e+09,'s^-1'), n=0.933333, Ea=(150.345,'kJ/mol'), T0=(1,'K'), comment="""Estimated using an average for rate rule [R3H_DS;Cd_rad_out_singleNd;XH_out]
Euclidian distance = 0
Multiplied by reaction path degeneracy 2.0
family: intra_H_migration"""),
)

reaction(
    label = 'reaction47',
    reactants = ['[CH]=[C]CC([CH2])=C[O](19973)'],
    products = ['[CH]=[C]C=C([CH2])[CH]O(21903)'],
    transitionState = 'TS47',
    kinetics = Arrhenius(A=(2.22e+06,'s^-1','*|/',3), n=1.78, Ea=(113.721,'kJ/mol'), T0=(1,'K'), Tmin=(600,'K'), Tmax=(2000,'K'), comment="""Estimated using an average for rate rule [R4H_SDS;O_rad_out;XH_out]
Euclidian distance = 0
Multiplied by reaction path degeneracy 2.0
family: intra_H_migration"""),
)

reaction(
    label = 'reaction48',
    reactants = ['[CH2][C]=CC([CH2])=[C]O(28971)'],
    products = ['[CH]=[C]C=C([CH2])[CH]O(21903)'],
    transitionState = 'TS48',
    kinetics = Arrhenius(A=(2.74e+08,'s^-1'), n=1.713, Ea=(181.895,'kJ/mol'), T0=(1,'K'), comment="""Estimated using template [RnH;Cd_rad_out_singleNd;Cd_H_out_singleH] for rate rule [R5Hall;Cd_rad_out_singleNd;Cd_H_out_singleH]
Euclidian distance = 1.0
Multiplied by reaction path degeneracy 2.0
family: intra_H_migration"""),
)

reaction(
    label = 'reaction49',
    reactants = ['[CH]=[C]C=C([CH2])[CH]O(21903)'],
    products = ['[CH2][C]1C=[C][CH]C1O(29849)'],
    transitionState = 'TS49',
    kinetics = Arrhenius(A=(9.13977e+10,'s^-1'), n=0.0396934, Ea=(32.3534,'kJ/mol'), T0=(1,'K'), comment="""Estimated using an average for rate rule [R5_linear;doublebond_intra;radadd_intra_cdsingleH]
Euclidian distance = 0
family: Intra_R_Add_Endocyclic"""),
)

reaction(
    label = 'reaction50',
    reactants = ['[CH]=[C]C=C([CH2])[CH]O(21903)'],
    products = ['[CH]=[C]C1CC1=CO(29817)'],
    transitionState = 'TS50',
    kinetics = Arrhenius(A=(5.94212e+13,'s^-1'), n=0.0123667, Ea=(5.39457,'kJ/mol'), T0=(1,'K'), comment="""Estimated using average of templates [Rn;Y_rad_out;Cpri_rad_out_2H] + [R3_SS;Y_rad_out;Ypri_rad_out] for rate rule [R3_SS;Y_rad_out;Cpri_rad_out_2H]
Euclidian distance = 2.0
family: Birad_recombination"""),
)

reaction(
    label = 'reaction51',
    reactants = ['[CH][C]([CH]C#C)CO(29850)'],
    products = ['[CH]=[C]C=C([CH2])[CH]O(21903)'],
    transitionState = 'TS51',
    kinetics = Arrhenius(A=(1.846e+10,'s^-1'), n=0.74, Ea=(145.185,'kJ/mol'), T0=(1,'K'), Tmin=(300,'K'), Tmax=(1500,'K'), comment="""Estimated using template [R3H_DS;Cd_rad_out_singleH;Cs_H_out_1H] for rate rule [R3H_DS;Cd_rad_out_singleH;Cs_H_out_H/NonDeO]
Euclidian distance = 1.0
Multiplied by reaction path degeneracy 2.0
family: intra_H_migration"""),
)

reaction(
    label = 'reaction52',
    reactants = ['[CH]C=CC([CH])=CO(29851)'],
    products = ['[CH]=[C]C=C([CH2])[CH]O(21903)'],
    transitionState = 'TS52',
    kinetics = Arrhenius(A=(1.13764e+12,'s^-1'), n=1.09983, Ea=(403.433,'kJ/mol'), T0=(1,'K'), comment="""Estimated using template [R4H_DSD;Cd_rad_out_single;XH_out] for rate rule [R4H_DSD;Cd_rad_out_singleH;XH_out]
Euclidian distance = 1.0
family: intra_H_migration"""),
)

reaction(
    label = 'reaction53',
    reactants = ['[CH]=[C]C=C([CH2])[CH]O(21903)'],
    products = ['[CH][C]1C=C([CH]O)C1(29852)'],
    transitionState = 'TS53',
    kinetics = Arrhenius(A=(1.953e+11,'s^-1'), n=0.387, Ea=(138.055,'kJ/mol'), T0=(1,'K'), comment="""Estimated using an average for rate rule [R4_D_D;doublebond_intra;radadd_intra]
Euclidian distance = 0
family: Intra_R_Add_Endocyclic"""),
)

reaction(
    label = 'reaction35',
    reactants = ['[CH]=[C][CH]C(O)[C]=C(22023)'],
    products = ['[CH]=[C]C=C([CH2])[CH]O(21903)'],
    transitionState = 'TS54',
    kinetics = Arrhenius(A=(8.66e+11,'s^-1'), n=0.438, Ea=(94.4747,'kJ/mol'), T0=(1,'K'), comment="""Estimated using an average for rate rule [cCs(-HR!H)CJ;CdsJ;C]
Euclidian distance = 0
family: 1,2_shiftC"""),
)

reaction(
    label = 'reaction55',
    reactants = ['[CH]C(=CO)C[C]=[CH](28754)'],
    products = ['[CH]=[C]C=C([CH2])[CH]O(21903)'],
    transitionState = 'TS55',
    kinetics = Arrhenius(A=(13437.7,'s^-1'), n=2.58467, Ea=(192.129,'kJ/mol'), T0=(1,'K'), comment="""Estimated using an average for rate rule [R3H_DS;Cd_rad_out_singleH;XH_out]
Euclidian distance = 0
Multiplied by reaction path degeneracy 2.0
family: intra_H_migration"""),
)

reaction(
    label = 'reaction56',
    reactants = ['[CH]=[C]C=C([CH2])[CH]O(21903)'],
    products = ['[CH]C(C=[C][CH2])=CO(24296)'],
    transitionState = 'TS56',
    kinetics = Arrhenius(A=(1.456e+11,'s^-1'), n=0.86, Ea=(191.209,'kJ/mol'), T0=(1,'K'), Tmin=(300,'K'), Tmax=(1500,'K'), comment="""Estimated using template [RnH;Cd_rad_out_singleH;Cd_H_out_singleH] for rate rule [R5Hall;Cd_rad_out_singleH;Cd_H_out_singleH]
Euclidian distance = 1.0
Multiplied by reaction path degeneracy 2.0
family: intra_H_migration"""),
)

reaction(
    label = 'reaction57',
    reactants = ['[CH]=[C]C=C([CH2])[CH]O(21903)'],
    products = ['O[CH][C]1[CH][C]=CC1(29853)'],
    transitionState = 'TS57',
    kinetics = Arrhenius(A=(9.13977e+10,'s^-1'), n=0.0396934, Ea=(32.3534,'kJ/mol'), T0=(1,'K'), comment="""Estimated using an average for rate rule [R5_linear;doublebond_intra;radadd_intra_cdsingleH]
Euclidian distance = 0
family: Intra_R_Add_Endocyclic"""),
)

reaction(
    label = 'reaction58',
    reactants = ['[CH]=[C]C=C([CH2])[CH]O(21903)'],
    products = ['[CH]=[C]CC(=C)C=O(19965)'],
    transitionState = 'TS58',
    kinetics = Arrhenius(A=(5.14222e+08,'s^-1'), n=0.311, Ea=(24.9733,'kJ/mol'), T0=(1,'K'), comment="""Estimated using template [R4;Y_rad;XH_Rrad] for rate rule [R4radEndo;Y_rad;XH_Rrad]
Euclidian distance = 1.0
family: Intra_Disproportionation"""),
)

reaction(
    label = 'reaction59',
    reactants = ['[CH]=[C]C=C([CH2])[CH]O(21903)'],
    products = ['[CH]=[C]C1C(=C)C1O(29824)'],
    transitionState = 'TS59',
    kinetics = Arrhenius(A=(1.61967e+11,'s^-1'), n=0.0247333, Ea=(23.7316,'kJ/mol'), T0=(1,'K'), comment="""Estimated using template [Rn;Y_rad_out;Cpri_rad_out_single] for rate rule [R3_SS;Y_rad_out;Cpri_rad_out_H/NonDeO]
Euclidian distance = 2.82842712475
family: Birad_recombination"""),
)

reaction(
    label = 'reaction60',
    reactants = ['[C]=[C]CC([CH2])=CO(28755)'],
    products = ['[CH]=[C]C=C([CH2])[CH]O(21903)'],
    transitionState = 'TS60',
    kinetics = Arrhenius(A=(2.30972e+07,'s^-1'), n=1.70216, Ea=(184.173,'kJ/mol'), T0=(1,'K'), comment="""Estimated using template [R3H_MS;Y_rad_out;XH_out] for rate rule [R3H_TS;Ct_rad_out;XH_out]
Euclidian distance = 1.41421356237
Multiplied by reaction path degeneracy 2.0
family: intra_H_migration"""),
)

reaction(
    label = 'reaction61',
    reactants = ['[C]#C[CH]C([CH2])[CH]O(29854)'],
    products = ['[CH]=[C]C=C([CH2])[CH]O(21903)'],
    transitionState = 'TS61',
    kinetics = Arrhenius(A=(322699,'s^-1'), n=1.96429, Ea=(152.474,'kJ/mol'), T0=(1,'K'), comment="""Estimated using template [R4Hall;Y_rad_out;XH_out] for rate rule [R4HJ_2;Ct_rad_out;XH_out]
Euclidian distance = 1.41421356237
family: intra_H_migration"""),
)

reaction(
    label = 'reaction62',
    reactants = ['[C]=[C]C=C([CH2])CO(29855)'],
    products = ['[CH]=[C]C=C([CH2])[CH]O(21903)'],
    transitionState = 'TS62',
    kinetics = Arrhenius(A=(26847.4,'s^-1'), n=2.09453, Ea=(64.7206,'kJ/mol'), T0=(1,'K'), comment="""Estimated using template [R5Hall;Y_rad_out;Cs_H_out_H/NonDeO] for rate rule [R5Hall;Ct_rad_out;Cs_H_out_H/NonDeO]
Euclidian distance = 1.0
Multiplied by reaction path degeneracy 2.0
family: intra_H_migration"""),
)

reaction(
    label = 'reaction63',
    reactants = ['[C][C]=CC(C)=CO(29856)'],
    products = ['[CH]=[C]C=C([CH2])[CH]O(21903)'],
    transitionState = 'TS63',
    kinetics = Arrhenius(A=(505536,'s^-1'), n=1.7378, Ea=(41.5716,'kJ/mol'), T0=(1,'K'), comment="""Estimated using template [R5Hall;Y_rad_out;Cs_H_out_2H] for rate rule [R5Hall;Ct_rad_out;Cs_H_out_2H]
Euclidian distance = 1.0
Multiplied by reaction path degeneracy 3.0
family: intra_H_migration"""),
)

reaction(
    label = 'reaction64',
    reactants = ['[C]#C(5143)', '[CH][C]([CH2])[CH]O(29857)'],
    products = ['[CH]=[C]C=C([CH2])[CH]O(21903)'],
    transitionState = 'TS64',
    kinetics = Arrhenius(A=(1.14854e+06,'m^3/(mol*s)'), n=0.575199, Ea=(34.3157,'kJ/mol'), T0=(1,'K'), comment="""Estimated using template [Y_rad;Birad] for rate rule [Ct_rad;Birad]
Euclidian distance = 1.0
family: Birad_R_Recombination"""),
)

reaction(
    label = 'reaction65',
    reactants = ['[CH]=C=[C]C([CH2])[CH]O(29858)'],
    products = ['[CH]=[C]C=C([CH2])[CH]O(21903)'],
    transitionState = 'TS65',
    kinetics = Arrhenius(A=(3.9024e+09,'s^-1'), n=1.09833, Ea=(177.686,'kJ/mol'), T0=(1,'K'), comment="""Estimated using an average for rate rule [R2H_S;Cd_rad_out_double;XH_out]
Euclidian distance = 0
family: intra_H_migration"""),
)

reaction(
    label = 'reaction66',
    reactants = ['[CH]=[C]C=C([CH2])[CH]O(21903)'],
    products = ['[CH2]C#CC([CH2])=CO(28959)'],
    transitionState = 'TS66',
    kinetics = Arrhenius(A=(7.437e+08,'s^-1'), n=1.045, Ea=(63.4002,'kJ/mol'), T0=(1,'K'), comment="""Estimated using an average for rate rule [R3radExo;Y_rad;XH_Rrad]
Euclidian distance = 0
family: Intra_Disproportionation"""),
)

reaction(
    label = 'reaction67',
    reactants = ['[CH]=C=[CH](18734)', '[CH2][C][CH]O(28660)'],
    products = ['[CH]=[C]C=C([CH2])[CH]O(21903)'],
    transitionState = 'TS67',
    kinetics = Arrhenius(A=(2.29708e+06,'m^3/(mol*s)'), n=0.575199, Ea=(34.3157,'kJ/mol'), T0=(1,'K'), comment="""Estimated using template [Y_rad;Birad] for rate rule [Cd_pri_rad;Birad]
Euclidian distance = 2.0
Multiplied by reaction path degeneracy 2.0
family: Birad_R_Recombination"""),
)

network(
    label = '5247',
    isomers = [
        '[CH]=[C]C=C([CH2])[CH]O(21903)',
    ],
    reactants = [
        ('C=C=CO(12571)', '[CH]=C=[CH](18734)'),
    ],
    bathGas = {
        'N2': 0.25,
        'Ne': 0.25,
        'He': 0.25,
        'Ar': 0.25,
    },
)

pressureDependence(
    label = '5247',
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

