species(
    label = '[CH2]C(=C)C=C=C[O](22603)',
    structure = SMILES('[CH2]C(=C)C=C=C[O]'),
    E0 = (281.078,'kJ/mol'),
    modes = [
        HarmonicOscillator(frequencies=([3000,3100,440,815,1455,1000,2950,3100,1380,975,1025,1650,540,610,2055,350,440,435,1725,2995,3025,975,1000,1300,1375,400,500,1630,1680,180,180],'cm^-1')),
        HinderedRotor(inertia=(1.43154,'amu*angstrom^2'), symmetry=1, barrier=(32.914,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(1.4384,'amu*angstrom^2'), symmetry=1, barrier=(33.0715,'kJ/mol'), semiclassical=False),
    ],
    spinMultiplicity = 3,
    opticalIsomers = 1,
    molecularWeight = (94.1112,'amu'),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[0.612281,0.0583592,-1.49787e-05,-3.75866e-08,2.29716e-11,33942.8,22.4721], Tmin=(100,'K'), Tmax=(929.873,'K')), NASAPolynomial(coeffs=[21.882,0.0082147,-7.92857e-07,6.55313e-11,-8.55619e-15,28199.5,-88.2081], Tmin=(929.873,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(281.078,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(299.321,'J/(mol*K)'), comment="""Thermo group additivity estimation: group(O2s-(Cds-Cd)H) + group(Cs-(Cds-Cds)HHH) + group(Cds-Cds(Cds-Cds)Cs) + group(Cds-Cds(Cds-Cds)H) + group(Cds-CdsHH) + group(Cds-CdsOsH) + group(Cdd-CdsCds) + radical(C=COJ) + radical(Allyl_P)"""),
)

species(
    label = 'C=C=C(6077)',
    structure = SMILES('C=C=C'),
    E0 = (182.728,'kJ/mol'),
    modes = [
        HarmonicOscillator(frequencies=([2950,3000,3050,3100,1330,1430,900,1050,1000,1050,1600,1700,540,610,2055],'cm^-1')),
    ],
    spinMultiplicity = 1,
    opticalIsomers = 1,
    molecularWeight = (40.0639,'amu'),
    collisionModel = TransportData(shapeIndex=2, epsilon=(2746.46,'J/mol'), sigma=(4.78521,'angstroms'), dipoleMoment=(0,'C*m'), polarizability=(0,'angstroms^3'), rotrelaxcollnum=0, comment="""Epsilon & sigma estimated with Tc=428.99 K, Pc=56.87 bar (from Joback method)"""),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[3.36163,0.00777114,2.52438e-05,-3.61894e-08,1.38594e-11,22005.7,7.1245], Tmin=(100,'K'), Tmax=(966.19,'K')), NASAPolynomial(coeffs=[6.46487,0.0106812,-3.73734e-06,6.87011e-10,-4.98581e-14,20670.6,-11.5463], Tmin=(966.19,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(182.728,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(157.975,'J/(mol*K)'), comment="""Thermo group additivity estimation: group(Cds-CdsHH) + group(Cds-CdsHH) + group(Cdd-CdsCds)"""),
)

species(
    label = 'C#CC=O(21959)',
    structure = SMILES('C#CC=O'),
    E0 = (84.2941,'kJ/mol'),
    modes = [
        HarmonicOscillator(frequencies=([2175,525,750,770,3400,2100,346.083,346.148,1254.41,1569.49,1569.49],'cm^-1')),
        HinderedRotor(inertia=(0.289847,'amu*angstrom^2'), symmetry=1, barrier=(24.6472,'kJ/mol'), semiclassical=False),
    ],
    spinMultiplicity = 1,
    opticalIsomers = 1,
    molecularWeight = (54.0474,'amu'),
    collisionModel = TransportData(shapeIndex=2, epsilon=(3153.34,'J/mol'), sigma=(5.09602,'angstroms'), dipoleMoment=(0,'C*m'), polarizability=(0,'angstroms^3'), rotrelaxcollnum=0, comment="""Epsilon & sigma estimated with Tc=492.54 K, Pc=54.07 bar (from Joback method)"""),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[3.00307,0.020875,-1.61807e-05,6.30357e-09,-1.00066e-12,10174.9,9.76778], Tmin=(100,'K'), Tmax=(1462.04,'K')), NASAPolynomial(coeffs=[7.33456,0.00902438,-4.02236e-06,7.59526e-10,-5.26541e-14,8908.32,-12.7744], Tmin=(1462.04,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(84.2941,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(128.874,'J/(mol*K)'), comment="""Thermo group additivity estimation: group(Cds-O2d(Cds-Cds)H) + group(Ct-CtCs) + group(Ct-CtH)"""),
)

species(
    label = 'C=C1CC1[C]=C[O](24299)',
    structure = SMILES('C=C1CC1[C]=C[O]'),
    E0 = (429.158,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
    molecularWeight = (94.1112,'amu'),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[1.05315,0.0528755,-1.80917e-05,-2.13437e-08,1.39419e-11,51732.7,21.4416], Tmin=(100,'K'), Tmax=(956.603,'K')), NASAPolynomial(coeffs=[16.7879,0.0157899,-4.95593e-06,8.74147e-10,-6.34439e-14,47408.8,-60.6368], Tmin=(956.603,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(429.158,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(303.478,'J/(mol*K)'), comment="""Thermo group additivity estimation: group(O2s-(Cds-Cd)H) + group(Cs-(Cds-Cds)(Cds-Cds)CsH) + group(Cs-(Cds-Cds)CsHH) + group(Cds-CdsCsCs) + group(Cds-CdsCsH) + group(Cds-CdsOsH) + group(Cds-CdsHH) + ring(Methylene_cyclopropane) + radical(C=COJ) + radical(Cds_S)"""),
)

species(
    label = '[CH2]C1([CH2])C=C=CO1(24365)',
    structure = SMILES('[CH2]C1([CH2])C=C=CO1'),
    E0 = (626.123,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
    molecularWeight = (94.1112,'amu'),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[-1.11237,0.082931,-8.72759e-05,4.24815e-08,-7.5302e-12,75514.8,25.5442], Tmin=(100,'K'), Tmax=(1655.98,'K')), NASAPolynomial(coeffs=[23.7138,0.00354968,2.21351e-06,-6.2475e-10,4.62265e-14,69954.4,-98.7127], Tmin=(1655.98,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(626.123,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(299.321,'J/(mol*K)'), comment="""Thermo group additivity estimation: group(O2s-Cs(Cds-Cd)) + group(Cs-(Cds-Cds)CsCsOs) + group(Cs-CsHHH) + group(Cs-CsHHH) + group(Cds-CdsCsH) + group(Cds-CdsOsH) + group(Cdd-CdsCds) + ring(1,2-Cyclopentadiene) + radical(CJC(C)OC) + radical(CJC(C)OC)"""),
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
    label = '[CH2]C(=C)C=C=C=O(24366)',
    structure = SMILES('[CH2]C(=C)C=C=C=O'),
    E0 = (325.175,'kJ/mol'),
    modes = [
        HarmonicOscillator(frequencies=([540,610,2055,350,440,435,1725,2950,3100,1380,975,1025,1650,3000,3100,440,815,1455,1000,3010,987.5,1337.5,450,1655,2120,512.5,787.5,180],'cm^-1')),
        HinderedRotor(inertia=(0.00039259,'amu*angstrom^2'), symmetry=1, barrier=(4.45749,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.000391138,'amu*angstrom^2'), symmetry=1, barrier=(4.44101,'kJ/mol'), semiclassical=False),
    ],
    spinMultiplicity = 2,
    opticalIsomers = 1,
    molecularWeight = (93.1033,'amu'),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[1.89123,0.0348035,3.30808e-06,-3.54943e-08,1.78267e-11,39196,5.08332], Tmin=(100,'K'), Tmax=(949.344,'K')), NASAPolynomial(coeffs=[14.677,0.00968133,-2.42374e-06,4.30569e-10,-3.42391e-14,35472.8,-62.7592], Tmin=(949.344,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(325.175,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(274.378,'J/(mol*K)'), comment="""Thermo group additivity estimation: group(Cs-(Cds-Cds)HHH) + group(Cds-Cds(Cds-Cds)Cs) + group(Cds-Cds(Cds-Cds)H) + group(Cds-CdsHH) + group(Cdd-CdsCds) + radical(Allyl_P)"""),
)

species(
    label = '[CH]=C=C[O](8556)',
    structure = SMILES('[CH]=C=C[O]'),
    E0 = (269.875,'kJ/mol'),
    modes = [
        HarmonicOscillator(frequencies=([3120,650,792.5,1650,3010,987.5,1337.5,450,1655,540,610,2055],'cm^-1')),
    ],
    spinMultiplicity = 3,
    opticalIsomers = 1,
    molecularWeight = (54.0474,'amu'),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[2.7401,0.0176211,1.62543e-05,-4.58526e-08,2.291e-11,32513.4,12.64], Tmin=(100,'K'), Tmax=(883.628,'K')), NASAPolynomial(coeffs=[13.5244,-0.00323064,4.17673e-06,-9.22668e-10,6.45049e-14,29515.7,-44.2318], Tmin=(883.628,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(269.875,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(133.032,'J/(mol*K)'), comment="""Thermo group additivity estimation: group(O2s-(Cds-Cd)H) + group(Cds-CdsOsH) + group(Cds-CdsHH) + group(Cdd-CdsCds) + radical(C=C=CJ) + radical(C=COJ)"""),
)

species(
    label = '[CH2]C([CH2])=CC#CO(24367)',
    structure = SMILES('[CH2]C([CH2])=CC#CO'),
    E0 = (296.994,'kJ/mol'),
    modes = [
        HarmonicOscillator(frequencies=([3000,3033.33,3066.67,3100,415,465,780,850,1435,1475,900,1100,3010,987.5,1337.5,450,1655,3615,1277.5,1000,2100,2250,500,550,350,440,435,1725,731.375],'cm^-1')),
        HinderedRotor(inertia=(0.0414153,'amu*angstrom^2'), symmetry=1, barrier=(86.6668,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.865946,'amu*angstrom^2'), symmetry=1, barrier=(19.9098,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(3.76129,'amu*angstrom^2'), symmetry=1, barrier=(86.4794,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(3.76577,'amu*angstrom^2'), symmetry=1, barrier=(86.5824,'kJ/mol'), semiclassical=False),
    ],
    spinMultiplicity = 3,
    opticalIsomers = 1,
    molecularWeight = (94.1112,'amu'),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[1.26632,0.0572783,-4.89995e-05,2.28759e-08,-4.39485e-12,35820.9,21.902], Tmin=(100,'K'), Tmax=(1233.93,'K')), NASAPolynomial(coeffs=[11.1765,0.0251526,-9.94667e-06,1.77656e-09,-1.20022e-13,33375.2,-27.9922], Tmin=(1233.93,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(296.994,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(291.007,'J/(mol*K)'), comment="""Thermo group additivity estimation: group(O2s-CtH) + group(Cs-(Cds-Cds)HHH) + group(Cs-(Cds-Cds)HHH) + group(Cds-CdsCsCs) + group(Cds-CdsCtH) + group(Ct-Ct(Cds-Cds)) + group(Ct-CtOs) + radical(CTCC=CCJ) + radical(CTCC=CCJ)"""),
)

species(
    label = '[CH2][C](C)C#CC=O(24368)',
    structure = SMILES('[CH2][C](C)C#CC=O'),
    E0 = (326.918,'kJ/mol'),
    modes = [
        HarmonicOscillator(frequencies=([2750,2800,2850,1350,1500,750,1050,1375,1000,2782.5,750,1395,475,1775,1000,3000,3100,440,815,1455,1000,2100,2250,500,550,360,370,350,180],'cm^-1')),
        HinderedRotor(inertia=(0.222174,'amu*angstrom^2'), symmetry=1, barrier=(5.10822,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.222158,'amu*angstrom^2'), symmetry=1, barrier=(5.10784,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.0726339,'amu*angstrom^2'), symmetry=1, barrier=(63.1849,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(2.74811,'amu*angstrom^2'), symmetry=1, barrier=(63.1845,'kJ/mol'), semiclassical=False),
    ],
    spinMultiplicity = 3,
    opticalIsomers = 1,
    molecularWeight = (94.1112,'amu'),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[1.19709,0.0620783,-6.94706e-05,4.65366e-08,-1.27088e-11,39419.6,26.6636], Tmin=(100,'K'), Tmax=(949.668,'K')), NASAPolynomial(coeffs=[8.93354,0.0268397,-9.62121e-06,1.581e-09,-9.99754e-14,38069.8,-9.631], Tmin=(949.668,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(326.918,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(291.007,'J/(mol*K)'), comment="""Thermo group additivity estimation: group(Cs-CtCsCsH) + group(Cs-CsHHH) + group(Cs-CsHHH) + group(Cds-O2d(Cds-Cds)H) + group(Ct-CtCs) + group(Ct-CtCs) + radical(Isobutyl) + radical(Tert_Propargyl)"""),
)

species(
    label = '[CH]=C(C)C=C=C[O](24369)',
    structure = SMILES('[CH]=C(C)C=C=C[O]'),
    E0 = (376.675,'kJ/mol'),
    modes = [
        HarmonicOscillator(frequencies=([3120,650,792.5,1650,2750,2800,2850,1350,1500,750,1050,1375,1000,540,610,2055,350,440,435,1725,2995,3025,975,1000,1300,1375,400,500,1630,1680,180],'cm^-1')),
        HinderedRotor(inertia=(1.20625,'amu*angstrom^2'), symmetry=1, barrier=(27.7341,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(1.20196,'amu*angstrom^2'), symmetry=1, barrier=(27.6354,'kJ/mol'), semiclassical=False),
    ],
    spinMultiplicity = 3,
    opticalIsomers = 1,
    molecularWeight = (94.1112,'amu'),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[0.532152,0.0644673,-4.0594e-05,-6.14576e-09,1.0651e-11,45439.3,23.575], Tmin=(100,'K'), Tmax=(937.856,'K')), NASAPolynomial(coeffs=[19.9069,0.0113033,-2.69837e-06,4.21631e-10,-3.09971e-14,40509.1,-75.5637], Tmin=(937.856,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(376.675,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(299.321,'J/(mol*K)'), comment="""Thermo group additivity estimation: group(O2s-(Cds-Cd)H) + group(Cs-(Cds-Cds)HHH) + group(Cds-Cds(Cds-Cds)Cs) + group(Cds-Cds(Cds-Cds)H) + group(Cds-CdsHH) + group(Cds-CdsOsH) + group(Cdd-CdsCds) + radical(C=COJ) + radical(Cds_P)"""),
)

species(
    label = '[CH2]C([CH2])=C=C=CO(24370)',
    structure = SMILES('[CH2]C([CH2])=C=C=CO'),
    E0 = (342.432,'kJ/mol'),
    modes = [
        HarmonicOscillator(frequencies=([3000,3033.33,3066.67,3100,415,465,780,850,1435,1475,900,1100,540,563.333,586.667,610,1970,2140,3010,987.5,1337.5,450,1655,3615,1277.5,1000,350,440,435,1725],'cm^-1')),
        HinderedRotor(inertia=(1.4533,'amu*angstrom^2'), symmetry=1, barrier=(33.4143,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(1.46057,'amu*angstrom^2'), symmetry=1, barrier=(33.5813,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(1.44657,'amu*angstrom^2'), symmetry=1, barrier=(33.2596,'kJ/mol'), semiclassical=False),
    ],
    spinMultiplicity = 3,
    opticalIsomers = 1,
    molecularWeight = (94.1112,'amu'),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[-0.844616,0.0816896,-8.48711e-05,4.12494e-08,-7.38882e-12,41380.8,27.4695], Tmin=(100,'K'), Tmax=(1603.89,'K')), NASAPolynomial(coeffs=[23.4832,0.00559329,7.21284e-07,-3.23462e-10,2.57724e-14,35560.9,-95.2068], Tmin=(1603.89,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(342.432,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(295.164,'J/(mol*K)'), comment="""Thermo group additivity estimation: group(O2s-(Cds-Cd)H) + group(Cs-(Cds-Cds)HHH) + group(Cs-(Cds-Cds)HHH) + group(Cds-CdsCsCs) + group(Cds-CdsOsH) + group(Cdd-CdsCds) + group(Cdd-CdsCds) + radical(Allyl_P) + radical(Allyl_P)"""),
)

species(
    label = '[CH2]C(C)=C[C]=C=O(24371)',
    structure = SMILES('[CH2]C(C)=C[C]=C=O'),
    E0 = (309.231,'kJ/mol'),
    modes = [
        HarmonicOscillator(frequencies=([2750,2800,2850,1350,1500,750,1050,1375,1000,350,440,435,1725,3000,3100,440,815,1455,1000,1685,370,3010,987.5,1337.5,450,1655,2120,512.5,787.5,180],'cm^-1')),
        HinderedRotor(inertia=(0.111266,'amu*angstrom^2'), symmetry=1, barrier=(2.55824,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.30279,'amu*angstrom^2'), symmetry=1, barrier=(6.96173,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(2.52781,'amu*angstrom^2'), symmetry=1, barrier=(58.1192,'kJ/mol'), semiclassical=False),
    ],
    spinMultiplicity = 3,
    opticalIsomers = 1,
    molecularWeight = (94.1112,'amu'),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[1.97505,0.052056,-2.7491e-05,-3.39056e-08,4.432e-11,37257.2,22.2891], Tmin=(100,'K'), Tmax=(501.717,'K')), NASAPolynomial(coeffs=[5.78192,0.0349161,-1.57443e-05,2.96836e-09,-2.05484e-13,36709,4.89181], Tmin=(501.717,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(309.231,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(295.164,'J/(mol*K)'), comment="""Thermo group additivity estimation: group(Cs-(Cds-Cds)HHH) + group(Cs-(Cds-Cds)HHH) + group(Cds-CdsCsCs) + group(Cd-Cd(CCO)H) + group(Cds-(Cdd-O2d)CsH) + radical(Cds_S) + radical(C=CC=CCJ)"""),
)

species(
    label = '[CH]C(=C)C=C=CO(24372)',
    structure = SMILES('[CH]C(=C)C=C=CO'),
    E0 = (358.801,'kJ/mol'),
    modes = [
        HarmonicOscillator(frequencies=([2950,3100,1380,975,1025,1650,540,610,2055,2995,3025,975,1000,1300,1375,400,500,1630,1680,3615,1277.5,1000,350,440,435,1725,180,180,180,180],'cm^-1')),
        HinderedRotor(inertia=(1.93553,'amu*angstrom^2'), symmetry=1, barrier=(44.5017,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(1.93404,'amu*angstrom^2'), symmetry=1, barrier=(44.4675,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(1.93481,'amu*angstrom^2'), symmetry=1, barrier=(44.4851,'kJ/mol'), semiclassical=False),
    ],
    spinMultiplicity = 3,
    opticalIsomers = 1,
    molecularWeight = (94.1112,'amu'),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[0.207434,0.06766,-2.61742e-05,-3.16235e-08,2.22853e-11,43305,24.0373], Tmin=(100,'K'), Tmax=(914.905,'K')), NASAPolynomial(coeffs=[22.7996,0.0106494,-1.17601e-06,5.42989e-11,-4.23318e-15,37423.2,-92.4999], Tmin=(914.905,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(358.801,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(295.164,'J/(mol*K)'), comment="""Thermo group additivity estimation: group(O2s-(Cds-Cd)H) + group(Cs-(Cds-Cds)HHH) + group(Cds-Cds(Cds-Cds)Cs) + group(Cds-Cds(Cds-Cds)H) + group(Cds-CdsHH) + group(Cds-CdsOsH) + group(Cdd-CdsCds) + radical(AllylJ2_triplet)"""),
)

species(
    label = '[CH2][C]=C(6078)',
    structure = SMILES('[CH2][C]=C'),
    E0 = (395.465,'kJ/mol'),
    modes = [
        HarmonicOscillator(frequencies=([2950,3100,1380,975,1025,1650,3000,3100,440,815,1455,1000,1685,370],'cm^-1')),
        HinderedRotor(inertia=(0.243012,'amu*angstrom^2'), symmetry=1, barrier=(30.4931,'kJ/mol'), semiclassical=False),
    ],
    spinMultiplicity = 3,
    opticalIsomers = 1,
    molecularWeight = (40.0639,'amu'),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[3.28772,0.0102496,1.79964e-05,-2.86352e-08,1.11956e-11,47593.9,10.4766], Tmin=(100,'K'), Tmax=(969.996,'K')), NASAPolynomial(coeffs=[6.37267,0.0109726,-3.91218e-06,7.11399e-10,-5.07602e-14,46362.9,-7.57277], Tmin=(969.996,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(395.465,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(153.818,'J/(mol*K)'), comment="""Thermo group additivity estimation: group(Cs-(Cds-Cds)HHH) + group(Cds-CdsCsH) + group(Cds-CdsHH) + radical(Allyl_P) + radical(Cds_S)"""),
)

species(
    label = '[CH2]C([CH2])=C=C=C[O](24373)',
    structure = SMILES('[CH2]C([CH2])=C=C=C[O]'),
    E0 = (483.894,'kJ/mol'),
    modes = [
        HarmonicOscillator(frequencies=([540,563.333,586.667,610,1970,2140,3010,987.5,1337.5,450,1655,350,440,435,1725,3000,3033.33,3066.67,3100,415,465,780,850,1435,1475,900,1100,180],'cm^-1')),
        HinderedRotor(inertia=(1.64465,'amu*angstrom^2'), symmetry=1, barrier=(37.8137,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(1.646,'amu*angstrom^2'), symmetry=1, barrier=(37.8447,'kJ/mol'), semiclassical=False),
    ],
    spinMultiplicity = 4,
    opticalIsomers = 1,
    molecularWeight = (93.1033,'amu'),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[0.625721,0.0625147,-3.99661e-05,-4.90277e-09,9.66178e-12,58331.2,23.5834], Tmin=(100,'K'), Tmax=(949.018,'K')), NASAPolynomial(coeffs=[19.6104,0.0107913,-2.93501e-06,5.00038e-10,-3.75203e-14,53453.7,-73.7262], Tmin=(949.018,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(483.894,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(274.378,'J/(mol*K)'), comment="""Thermo group additivity estimation: group(O2s-(Cds-Cd)H) + group(Cs-(Cds-Cds)HHH) + group(Cs-(Cds-Cds)HHH) + group(Cds-CdsCsCs) + group(Cds-CdsOsH) + group(Cdd-CdsCds) + group(Cdd-CdsCds) + radical(Allyl_P) + radical(C=COJ) + radical(Allyl_P)"""),
)

species(
    label = '[CH]C(=C)C=C=C[O](22618)',
    structure = SMILES('[CH]C(=C)C=C=C[O]'),
    E0 = (500.264,'kJ/mol'),
    modes = [
        HarmonicOscillator(frequencies=([2950,3100,1380,975,1025,1650,540,610,2055,350,440,435,1725,2995,3025,975,1000,1300,1375,400,500,1630,1680,180,180,180,180,180],'cm^-1')),
        HinderedRotor(inertia=(2.00191,'amu*angstrom^2'), symmetry=1, barrier=(46.0278,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(2.00636,'amu*angstrom^2'), symmetry=1, barrier=(46.1302,'kJ/mol'), semiclassical=False),
    ],
    spinMultiplicity = 4,
    opticalIsomers = 1,
    molecularWeight = (93.1033,'amu'),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[0.61675,0.0607309,-2.30026e-05,-2.49436e-08,1.72533e-11,60302.2,23.9773], Tmin=(100,'K'), Tmax=(936.072,'K')), NASAPolynomial(coeffs=[19.4917,0.0144273,-3.85094e-06,6.21069e-10,-4.4817e-14,55263.5,-73.8753], Tmin=(936.072,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(500.264,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(274.378,'J/(mol*K)'), comment="""Thermo group additivity estimation: group(O2s-(Cds-Cd)H) + group(Cs-(Cds-Cds)HHH) + group(Cds-Cds(Cds-Cds)Cs) + group(Cds-Cds(Cds-Cds)H) + group(Cds-CdsHH) + group(Cds-CdsOsH) + group(Cdd-CdsCds) + radical(C=COJ) + radical(AllylJ2_triplet)"""),
)

species(
    label = '[CH2]C([CH2])=C[C]=C=O(24374)',
    structure = SMILES('[CH2]C([CH2])=C[C]=C=O'),
    E0 = (427.287,'kJ/mol'),
    modes = [
        HarmonicOscillator(frequencies=([3000,3033.33,3066.67,3100,415,465,780,850,1435,1475,900,1100,1685,370,2120,512.5,787.5,350,440,435,1725,3010,987.5,1337.5,450,1655,180],'cm^-1')),
        HinderedRotor(inertia=(2.70308,'amu*angstrom^2'), symmetry=1, barrier=(62.1492,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.0402159,'amu*angstrom^2'), symmetry=1, barrier=(62.1531,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.040229,'amu*angstrom^2'), symmetry=1, barrier=(62.1064,'kJ/mol'), semiclassical=False),
    ],
    spinMultiplicity = 4,
    opticalIsomers = 1,
    molecularWeight = (93.1033,'amu'),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[1.79123,0.050661,-4.40726e-05,2.36769e-08,-5.62769e-12,51468.5,24.3975], Tmin=(100,'K'), Tmax=(973.665,'K')), NASAPolynomial(coeffs=[6.70048,0.0304928,-1.3002e-05,2.40276e-09,-1.6529e-13,50512.5,0.844248], Tmin=(973.665,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(427.287,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(270.22,'J/(mol*K)'), comment="""Thermo group additivity estimation: group(Cs-(Cds-Cds)HHH) + group(Cs-(Cds-Cds)HHH) + group(Cds-CdsCsCs) + group(Cd-Cd(CCO)H) + group(Cds-(Cdd-O2d)CsH) + radical(C=CC=CCJ) + radical(Cds_S) + radical(C=CC=CCJ)"""),
)

species(
    label = '[O]C=C=C[C]1CC1(24375)',
    structure = SMILES('[O]C=C=C[C]1CC1'),
    E0 = (319.552,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
    molecularWeight = (94.1112,'amu'),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[1.64891,0.0301433,5.60373e-05,-1.04109e-07,4.48911e-11,38537.7,20.4522], Tmin=(100,'K'), Tmax=(932.869,'K')), NASAPolynomial(coeffs=[19.0151,0.0107561,-1.3489e-06,1.89642e-10,-2.07037e-14,32901.2,-74.9675], Tmin=(932.869,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(319.552,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(303.478,'J/(mol*K)'), comment="""Thermo group additivity estimation: group(O2s-(Cds-Cd)H) + group(Cs-(Cds-Cdd-Cd)CsCsH) + group(Cs-CsCsHH) + group(Cs-CsCsHH) + group(Cds-CdsCsH) + group(Cds-CdsOsH) + group(Cdd-CdsCds) + ring(Cyclopropane) + radical(C=COJ) + radical(Allyl_T)"""),
)

species(
    label = '[CH2]C(=C)C=C1[CH]O1(24376)',
    structure = SMILES('[CH2]C(=C)C=C1[CH]O1'),
    E0 = (307.261,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
    molecularWeight = (94.1112,'amu'),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[0.99735,0.0375789,6.08874e-05,-1.28432e-07,5.8314e-11,37089.5,20.3649], Tmin=(100,'K'), Tmax=(920.643,'K')), NASAPolynomial(coeffs=[27.1708,-0.00225724,5.41659e-06,-1.09606e-09,6.56199e-14,29139.1,-120.748], Tmin=(920.643,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(307.261,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(299.321,'J/(mol*K)'), comment="""Thermo group additivity estimation: group(O2s-Cs(Cds-Cd)) + group(Cs-(Cds-Cds)OsHH) + group(Cs-(Cds-Cds)HHH) + group(Cds-Cds(Cds-Cds)Cs) + group(Cds-CdsCsOs) + group(Cds-Cds(Cds-Cds)H) + group(Cds-CdsHH) + ring(methyleneoxirane) + radical(Allyl_P) + radical(C=CCJO)"""),
)

species(
    label = 'C=C1[CH]C(=C[O])C1(24377)',
    structure = SMILES('C=C1[CH]C(=C[O])C1'),
    E0 = (230.531,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
    molecularWeight = (94.1112,'amu'),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[2.03083,0.0107456,0.000127144,-1.87761e-07,7.70381e-11,27827.5,21.7244], Tmin=(100,'K'), Tmax=(925.144,'K')), NASAPolynomial(coeffs=[23.1423,0.00325848,3.42628e-06,-7.08863e-10,3.61638e-14,20335.5,-97.8637], Tmin=(925.144,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(230.531,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(307.635,'J/(mol*K)'), comment="""Thermo group additivity estimation: group(O2s-(Cds-Cd)H) + group(Cs-(Cds-Cds)(Cds-Cds)HH) + group(Cs-(Cds-Cds)(Cds-Cds)HH) + group(Cds-CdsCsCs) + group(Cds-CdsCsCs) + group(Cds-CdsOsH) + group(Cds-CdsHH) + ring(Cyclobutane) + radical(C=CCJC=C) + radical(C=COJ)"""),
)

species(
    label = '[CH2]C(=C)C1[C]=CO1(24378)',
    structure = SMILES('[CH2]C(=C)C1[C]=CO1'),
    E0 = (413.034,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
    molecularWeight = (94.1112,'amu'),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[1.05157,0.0408095,4.10421e-05,-9.86001e-08,4.49665e-11,49804.8,21.7334], Tmin=(100,'K'), Tmax=(935.772,'K')), NASAPolynomial(coeffs=[23.858,0.00453805,1.05787e-06,-2.07273e-10,3.62719e-15,42856.3,-101.101], Tmin=(935.772,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(413.034,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(299.321,'J/(mol*K)'), comment="""Thermo group additivity estimation: group(O2s-Cs(Cds-Cd)) + group(Cs-(Cds-Cds)(Cds-Cds)OsH) + group(Cs-(Cds-Cds)HHH) + group(Cds-CdsCsCs) + group(Cds-CdsCsH) + group(Cds-CdsOsH) + group(Cds-CdsHH) + ring(Cyclobutene) + radical(Allyl_P) + radical(Cds_S)"""),
)

species(
    label = 'C=C1C=[C]C([O])C1(24348)',
    structure = SMILES('C=C1C=[C]C([O])C1'),
    E0 = (412.017,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
    molecularWeight = (94.1112,'amu'),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[1.7623,0.0323228,3.61748e-05,-7.24354e-08,3.03927e-11,49650,21.4299], Tmin=(100,'K'), Tmax=(968.401,'K')), NASAPolynomial(coeffs=[15.3207,0.0180053,-6.21704e-06,1.19847e-09,-9.16519e-14,45069.4,-53.6383], Tmin=(968.401,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(412.017,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(307.635,'J/(mol*K)'), comment="""Thermo group additivity estimation: group(O2s-CsH) + group(Cs-(Cds-Cds)CsOsH) + group(Cs-(Cds-Cds)CsHH) + group(Cds-Cds(Cds-Cds)Cs) + group(Cds-CdsCsH) + group(Cds-Cds(Cds-Cds)H) + group(Cds-CdsHH) + ring(3-Methylenecyclopentene) + radical(CC(C)OJ) + radical(cyclopentene-vinyl)"""),
)

species(
    label = '[CH2]C1=C[C]=COC1(24334)',
    structure = SMILES('[CH2]C1=C[C]=COC1'),
    E0 = (226.721,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
    molecularWeight = (94.1112,'amu'),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[1.56953,0.0309304,6.08583e-05,-1.13189e-07,4.94519e-11,27376.6,17.3443], Tmin=(100,'K'), Tmax=(919.94,'K')), NASAPolynomial(coeffs=[19.7886,0.0100141,-1.00658e-07,-1.21535e-10,3.31453e-15,21557.5,-82.4402], Tmin=(919.94,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(226.721,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(303.478,'J/(mol*K)'), comment="""Thermo group additivity estimation: group(O2s-Cs(Cds-Cd)) + group(Cs-(Cds-Cds)OsHH) + group(Cs-(Cds-Cds)HHH) + group(Cds-CdsCsCs) + group(Cds-Cds(Cds-Cds)H) + group(Cds-Cds(Cds-Cds)H) + group(Cds-CdsOsH) + ring(1,3-Cyclohexadiene) + radical(C=CJC=C) + radical(C=CC=CCJ)"""),
)

species(
    label = 'C=C(C)C=C=C=O(24379)',
    structure = SMILES('C=C(C)C=C=C=O'),
    E0 = (173.676,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
    molecularWeight = (94.1112,'amu'),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[1.87689,0.0374752,-5.90489e-06,-2.1193e-08,1.13667e-11,20973.1,5.54733], Tmin=(100,'K'), Tmax=(977.993,'K')), NASAPolynomial(coeffs=[12.5345,0.015481,-5.29397e-06,9.69218e-10,-7.02118e-14,17855.7,-50.9127], Tmin=(977.993,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(173.676,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(299.321,'J/(mol*K)'), comment="""Thermo group additivity estimation: group(Cs-(Cds-Cds)HHH) + group(Cds-Cds(Cds-Cds)Cs) + group(Cds-Cds(Cds-Cds)H) + group(Cds-CdsHH) + group(Cdd-CdsCds)"""),
)

species(
    label = '[CH2]C([CH2])[C]=C=C[O](23592)',
    structure = SMILES('[CH2]C([CH2])[C]=C=C[O]'),
    E0 = (672.981,'kJ/mol'),
    modes = [
        HarmonicOscillator(frequencies=([1685,370,540,610,2055,3010,987.5,1337.5,450,1655,1380,1390,370,380,2900,435,3000,3033.33,3066.67,3100,415,465,780,850,1435,1475,900,1100,368.949,372.833],'cm^-1')),
        HinderedRotor(inertia=(0.113988,'amu*angstrom^2'), symmetry=1, barrier=(11.1672,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(1.22919,'amu*angstrom^2'), symmetry=1, barrier=(119.627,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(1.19802,'amu*angstrom^2'), symmetry=1, barrier=(119.627,'kJ/mol'), semiclassical=False),
    ],
    spinMultiplicity = 5,
    opticalIsomers = 1,
    molecularWeight = (94.1112,'amu'),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[0.280306,0.0702613,-7.32263e-05,3.88442e-08,-7.83627e-12,81084.8,30.8976], Tmin=(100,'K'), Tmax=(1368.07,'K')), NASAPolynomial(coeffs=[16.8728,0.0138306,-2.67298e-06,2.33087e-10,-7.48051e-15,77285.7,-51.6436], Tmin=(1368.07,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(672.981,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(295.164,'J/(mol*K)'), comment="""Thermo group additivity estimation: group(O2s-(Cds-Cd)H) + group(Cs-(Cds-Cdd-Cd)CsCsH) + group(Cs-CsHHH) + group(Cs-CsHHH) + group(Cds-CdsCsH) + group(Cds-CdsOsH) + group(Cdd-CdsCds) + radical(Cds_S) + radical(Isobutyl) + radical(Isobutyl) + radical(C=COJ)"""),
)

species(
    label = '[CH2]C([CH2])=[C]C=C[O](24380)',
    structure = SMILES('[CH2]C([CH2])=[C]C=C[O]'),
    E0 = (386.618,'kJ/mol'),
    modes = [
        HarmonicOscillator(frequencies=([3000,3033.33,3066.67,3100,415,465,780,850,1435,1475,900,1100,1685,370,350,440,435,1725,2995,3025,975,1000,1300,1375,400,500,1630,1680,210.868,212.118],'cm^-1')),
        HinderedRotor(inertia=(1.54566,'amu*angstrom^2'), symmetry=1, barrier=(49.5448,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(1.56162,'amu*angstrom^2'), symmetry=1, barrier=(49.5072,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(1.54392,'amu*angstrom^2'), symmetry=1, barrier=(49.5104,'kJ/mol'), semiclassical=False),
    ],
    spinMultiplicity = 5,
    opticalIsomers = 1,
    molecularWeight = (94.1112,'amu'),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[1.05502,0.0499206,3.21042e-06,-5.2387e-08,2.80458e-11,46619.3,25.1125], Tmin=(100,'K'), Tmax=(903.85,'K')), NASAPolynomial(coeffs=[18.4768,0.0128069,-1.55775e-06,7.66653e-11,-3.76516e-15,41836.6,-66.2121], Tmin=(903.85,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(386.618,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(295.164,'J/(mol*K)'), comment="""Thermo group additivity estimation: group(O2s-(Cds-Cd)H) + group(Cs-(Cds-Cds)HHH) + group(Cs-(Cds-Cds)HHH) + group(Cds-CdsCsCs) + group(Cds-Cds(Cds-Cds)H) + group(Cds-Cds(Cds-Cds)H) + group(Cds-CdsOsH) + radical(C=CJC=C) + radical(C=CC=CCJ) + radical(C=CC=CCJ) + radical(C=COJ)"""),
)

species(
    label = '[CH2]C([CH2])=CC=[C][O](24381)',
    structure = SMILES('[CH2]C([CH2])=CC=[C][O]'),
    E0 = (427.366,'kJ/mol'),
    modes = [
        HarmonicOscillator(frequencies=([3000,3033.33,3066.67,3100,415,465,780,850,1435,1475,900,1100,1685,370,350,440,435,1725,2995,3025,975,1000,1300,1375,400,500,1630,1680,523.806,524.198],'cm^-1')),
        HinderedRotor(inertia=(0.275755,'amu*angstrom^2'), symmetry=1, barrier=(53.8342,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.276273,'amu*angstrom^2'), symmetry=1, barrier=(53.8386,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.2761,'amu*angstrom^2'), symmetry=1, barrier=(53.8266,'kJ/mol'), semiclassical=False),
    ],
    spinMultiplicity = 5,
    opticalIsomers = 1,
    molecularWeight = (94.1112,'amu'),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[1.33706,0.0452312,4.93061e-06,-4.55642e-08,2.30497e-11,51508.5,27.6304], Tmin=(100,'K'), Tmax=(928.599,'K')), NASAPolynomial(coeffs=[15.8191,0.0169404,-4.43935e-06,6.98423e-10,-4.91141e-14,47349,-49.0786], Tmin=(928.599,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(427.366,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(295.164,'J/(mol*K)'), comment="""Thermo group additivity estimation: group(O2s-(Cds-Cd)H) + group(Cs-(Cds-Cds)HHH) + group(Cs-(Cds-Cds)HHH) + group(Cds-CdsCsCs) + group(Cds-Cds(Cds-Cds)H) + group(Cds-Cds(Cds-Cds)H) + group(Cds-CdsOsH) + radical(C=CC=CCJ) + radical(C=CJO) + radical(C=COJ) + radical(C=CC=CCJ)"""),
)

species(
    label = '[CH]C(=C)C[C]=C[O](20209)',
    structure = SMILES('[CH]C(=C)C[C]=C[O]'),
    E0 = (590.891,'kJ/mol'),
    modes = [
        HarmonicOscillator(frequencies=([1685,370,2750,2850,1437.5,1250,1305,750,350,3010,987.5,1337.5,450,1655,350,440,435,1725,2950,3100,1380,975,1025,1650,276.118,277.905,278.594,278.905,278.908,280.512],'cm^-1')),
        HinderedRotor(inertia=(0.889306,'amu*angstrom^2'), symmetry=1, barrier=(49.3447,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.889516,'amu*angstrom^2'), symmetry=1, barrier=(49.3395,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.901456,'amu*angstrom^2'), symmetry=1, barrier=(49.3076,'kJ/mol'), semiclassical=False),
    ],
    spinMultiplicity = 5,
    opticalIsomers = 1,
    molecularWeight = (94.1112,'amu'),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[0.73953,0.0613204,-3.19944e-05,-6.15398e-09,7.86442e-12,71194.4,28.1347], Tmin=(100,'K'), Tmax=(982.344,'K')), NASAPolynomial(coeffs=[15.7452,0.0230487,-8.41514e-06,1.50184e-09,-1.047e-13,67144.7,-49.5984], Tmin=(982.344,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(590.891,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(295.164,'J/(mol*K)'), comment="""Thermo group additivity estimation: group(O2s-(Cds-Cd)H) + group(Cs-(Cds-Cds)(Cds-Cds)HH) + group(Cs-(Cds-Cds)HHH) + group(Cds-CdsCsCs) + group(Cds-CdsCsH) + group(Cds-CdsOsH) + group(Cds-CdsHH) + radical(Cds_S) + radical(AllylJ2_triplet) + radical(C=COJ)"""),
)

species(
    label = '[CH]C([CH2])=CC=C[O](20208)',
    structure = SMILES('[CH]C([CH2])=CC=C[O]'),
    E0 = (440.251,'kJ/mol'),
    modes = [
        HarmonicOscillator(frequencies=([3000,3100,440,815,1455,1000,350,440,435,1725,2995,3010,3025,975,987.5,1000,1300,1337.5,1375,400,450,500,1630,1655,1680,180,180,180,180,180],'cm^-1')),
        HinderedRotor(inertia=(2.13029,'amu*angstrom^2'), symmetry=1, barrier=(48.9796,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(2.13018,'amu*angstrom^2'), symmetry=1, barrier=(48.977,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(2.12911,'amu*angstrom^2'), symmetry=1, barrier=(48.9524,'kJ/mol'), semiclassical=False),
    ],
    spinMultiplicity = 5,
    opticalIsomers = 1,
    molecularWeight = (94.1112,'amu'),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[0.852906,0.0522324,9.17687e-06,-5.92609e-08,2.96226e-11,53078.9,25.7903], Tmin=(100,'K'), Tmax=(928.715,'K')), NASAPolynomial(coeffs=[19.0791,0.0171229,-4.19955e-06,6.49426e-10,-4.70293e-14,47822.2,-70.8674], Tmin=(928.715,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(440.251,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(295.164,'J/(mol*K)'), comment="""Thermo group additivity estimation: group(O2s-(Cds-Cd)H) + group(Cs-(Cds-Cds)HHH) + group(Cs-(Cds-Cds)HHH) + group(Cds-CdsCsCs) + group(Cds-Cds(Cds-Cds)H) + group(Cds-Cds(Cds-Cds)H) + group(Cds-CdsOsH) + radical(C=COJ) + radical(AllylJ2_triplet) + radical(C=CC=CCJ)"""),
)

species(
    label = '[CH2]C([CH2])[CH][C]=C=O(23593)',
    structure = SMILES('[CH2]C([CH2])[CH][C]=C=O'),
    E0 = (601.462,'kJ/mol'),
    modes = [
        HarmonicOscillator(frequencies=([1685,370,3025,407.5,1350,352.5,2120,512.5,787.5,1380,1390,370,380,2900,435,3000,3033.33,3066.67,3100,415,465,780,850,1435,1475,900,1100,353.873,354.029],'cm^-1')),
        HinderedRotor(inertia=(0.00134601,'amu*angstrom^2'), symmetry=1, barrier=(0.119627,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.00134535,'amu*angstrom^2'), symmetry=1, barrier=(0.119627,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.00134613,'amu*angstrom^2'), symmetry=1, barrier=(0.119627,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.625164,'amu*angstrom^2'), symmetry=1, barrier=(55.5657,'kJ/mol'), semiclassical=False),
    ],
    spinMultiplicity = 5,
    opticalIsomers = 1,
    molecularWeight = (94.1112,'amu'),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[1.12979,0.0657864,-7.198e-05,4.16895e-08,-8.22884e-12,72440.3,28.6589], Tmin=(100,'K'), Tmax=(730.007,'K')), NASAPolynomial(coeffs=[9.71738,0.027342,-1.06778e-05,1.86362e-09,-1.23276e-13,70957,-11.6402], Tmin=(730.007,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(601.462,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(291.007,'J/(mol*K)'), comment="""Thermo group additivity estimation: group(Cs-CsCsCsH) + group(Cs-(Cds-Cdd-O2d)CsHH) + group(Cs-CsHHH) + group(Cs-CsHHH) + group(Cds-(Cdd-O2d)CsH) + radical(Isobutyl) + radical(CCJC(C)=C=O) + radical(Isobutyl) + radical(CCCJ=C=O)"""),
)

species(
    label = '[CH]C(=C)C=[C]C[O](24382)',
    structure = SMILES('[CH]C(=C)C=[C]C[O]'),
    E0 = (701.401,'kJ/mol'),
    modes = [
        HarmonicOscillator(frequencies=([1685,370,2750,2850,1437.5,1250,1305,750,350,2950,3100,1380,975,1025,1650,350,440,435,1725,3010,987.5,1337.5,450,1655,428.285,428.339,428.373,428.417,428.423,428.46],'cm^-1')),
        HinderedRotor(inertia=(0.407331,'amu*angstrom^2'), symmetry=1, barrier=(53.0663,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.407579,'amu*angstrom^2'), symmetry=1, barrier=(53.065,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.407573,'amu*angstrom^2'), symmetry=1, barrier=(53.0667,'kJ/mol'), semiclassical=False),
    ],
    spinMultiplicity = 5,
    opticalIsomers = 1,
    molecularWeight = (94.1112,'amu'),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[1.5478,0.0607301,-4.4784e-05,7.90487e-10,1.63647e-11,84440.8,25.4766], Tmin=(100,'K'), Tmax=(557.075,'K')), NASAPolynomial(coeffs=[5.87991,0.0402715,-1.83669e-05,3.48692e-09,-2.43046e-13,83793,5.62833], Tmin=(557.075,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(701.401,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(295.164,'J/(mol*K)'), comment="""Thermo group additivity estimation: group(O2s-CsH) + group(Cs-(Cds-Cds)OsHH) + group(Cs-(Cds-Cds)HHH) + group(Cds-Cds(Cds-Cds)Cs) + group(Cds-CdsCsH) + group(Cds-Cds(Cds-Cds)H) + group(Cds-CdsHH) + radical(CCOJ) + radical(AllylJ2_triplet) + radical(Cds_S)"""),
)

species(
    label = '[CH2]C(C)=C[C]=[C][O](24383)',
    structure = SMILES('[CH2]C(C)=C[C]=[C][O]'),
    E0 = (508.306,'kJ/mol'),
    modes = [
        HarmonicOscillator(frequencies=([3000,3100,440,815,1455,1000,2750,2800,2850,1350,1500,750,1050,1375,1000,350,440,435,1725,1670,1700,300,440,3010,987.5,1337.5,450,1655,219.576,219.576],'cm^-1')),
        HinderedRotor(inertia=(2.28997,'amu*angstrom^2'), symmetry=1, barrier=(78.3476,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(2.28997,'amu*angstrom^2'), symmetry=1, barrier=(78.3476,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(2.28998,'amu*angstrom^2'), symmetry=1, barrier=(78.3475,'kJ/mol'), semiclassical=False),
    ],
    spinMultiplicity = 5,
    opticalIsomers = 1,
    molecularWeight = (94.1112,'amu'),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[0.665506,0.0643618,-6.02595e-05,2.98905e-08,-5.8394e-12,61262.8,28.0799], Tmin=(100,'K'), Tmax=(1308.37,'K')), NASAPolynomial(coeffs=[14.7826,0.0190331,-5.80458e-06,8.76216e-10,-5.32748e-14,57754.4,-43.1118], Tmin=(1308.37,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(508.306,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(295.164,'J/(mol*K)'), comment="""Thermo group additivity estimation: group(O2s-(Cds-Cd)H) + group(Cs-(Cds-Cds)HHH) + group(Cs-(Cds-Cds)HHH) + group(Cds-CdsCsCs) + group(Cds-Cds(Cds-Cds)H) + group(Cds-Cds(Cds-Cds)H) + group(Cds-CdsOsH) + radical(C=CJC=C) + radical(C=CJO) + radical(C=COJ) + radical(C=CC=CCJ)"""),
)

species(
    label = 'C=[C]CC=C=C[O](22607)',
    structure = SMILES('C=[C]CC=C=C[O]'),
    E0 = (399.84,'kJ/mol'),
    modes = [
        HarmonicOscillator(frequencies=([1685,370,2750,2850,1437.5,1250,1305,750,350,540,610,2055,2995,3025,975,1000,1300,1375,400,500,1630,1680,2950,3100,1380,975,1025,1650,180,180,180],'cm^-1')),
        HinderedRotor(inertia=(1.03654,'amu*angstrom^2'), symmetry=1, barrier=(23.832,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(1.03759,'amu*angstrom^2'), symmetry=1, barrier=(23.8561,'kJ/mol'), semiclassical=False),
    ],
    spinMultiplicity = 3,
    opticalIsomers = 1,
    molecularWeight = (94.1112,'amu'),
    collisionModel = TransportData(shapeIndex=2, epsilon=(3942.94,'J/mol'), sigma=(6.38069,'angstroms'), dipoleMoment=(0,'C*m'), polarizability=(0,'angstroms^3'), rotrelaxcollnum=0, comment="""Epsilon & sigma estimated with Tc=615.88 K, Pc=34.44 bar (from Joback method)"""),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[0.951029,0.0569117,-3.25147e-05,-4.32801e-09,7.28063e-12,48208.5,26.3697], Tmin=(100,'K'), Tmax=(980.593,'K')), NASAPolynomial(coeffs=[16.2574,0.0170218,-5.98614e-06,1.08493e-09,-7.75848e-14,44122.6,-52.7022], Tmin=(980.593,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(399.84,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(299.321,'J/(mol*K)'), comment="""Thermo group additivity estimation: group(O2s-(Cds-Cd)H) + group(Cs-(Cds-Cds)(Cds-Cds)HH) + group(Cds-CdsCsH) + group(Cds-CdsCsH) + group(Cds-CdsHH) + group(Cds-CdsOsH) + group(Cdd-CdsCds) + radical(Cds_S) + radical(C=COJ)"""),
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
    label = 'C=[C]C=C=C[O](23378)',
    structure = SMILES('C=[C]C=C=C[O]'),
    E0 = (366.165,'kJ/mol'),
    modes = [
        HarmonicOscillator(frequencies=([2950,3100,1380,975,1025,1650,540,610,2055,2995,3025,975,1000,1300,1375,400,500,1630,1680,1685,370,180,180],'cm^-1')),
        HinderedRotor(inertia=(1.63597,'amu*angstrom^2'), symmetry=1, barrier=(37.6142,'kJ/mol'), semiclassical=False),
    ],
    spinMultiplicity = 3,
    opticalIsomers = 1,
    molecularWeight = (80.0847,'amu'),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[0.442946,0.0581895,-5.88988e-05,2.83071e-08,-4.99578e-12,44184.9,22.5047], Tmin=(100,'K'), Tmax=(1656.64,'K')), NASAPolynomial(coeffs=[16.9916,0.00540583,5.07716e-07,-2.72871e-10,2.24426e-14,40461.9,-60.3744], Tmin=(1656.64,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(366.165,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(228.648,'J/(mol*K)'), comment="""Thermo group additivity estimation: group(O2s-(Cds-Cd)H) + group(Cds-Cds(Cds-Cds)H) + group(Cds-Cds(Cds-Cds)H) + group(Cds-CdsOsH) + group(Cds-CdsHH) + group(Cdd-CdsCds) + radical(C=COJ) + radical(C=CJC=C)"""),
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
    label = 'C#CC=C([CH2])[CH2](19657)',
    structure = SMILES('C#CC=C([CH2])[CH2]'),
    E0 = (438.46,'kJ/mol'),
    modes = [
        HarmonicOscillator(frequencies=([3000,3033.33,3066.67,3100,415,465,780,850,1435,1475,900,1100,750,770,3400,2100,350,440,435,1725,2175,525,3010,987.5,1337.5,450,1655],'cm^-1')),
        HinderedRotor(inertia=(0.142552,'amu*angstrom^2'), symmetry=1, barrier=(71.2344,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(3.10175,'amu*angstrom^2'), symmetry=1, barrier=(71.3153,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(3.10514,'amu*angstrom^2'), symmetry=1, barrier=(71.3934,'kJ/mol'), semiclassical=False),
    ],
    spinMultiplicity = 3,
    opticalIsomers = 1,
    molecularWeight = (78.1118,'amu'),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[1.39256,0.0509418,-4.21539e-05,1.91619e-08,-3.53282e-12,52833.8,18.1133], Tmin=(100,'K'), Tmax=(1303.18,'K')), NASAPolynomial(coeffs=[11.377,0.0202952,-6.87877e-06,1.11622e-09,-7.09586e-14,50231.5,-32.7], Tmin=(1303.18,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(438.46,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(270.22,'J/(mol*K)'), comment="""Thermo group additivity estimation: group(Cs-(Cds-Cds)HHH) + group(Cs-(Cds-Cds)HHH) + group(Cds-CdsCsCs) + group(Cds-CdsCtH) + group(Ct-Ct(Cds-Cds)) + group(Ct-CtH) + radical(CTCC=CCJ) + radical(CTCC=CCJ)"""),
)

species(
    label = '[CH2]C1([CH2])C=C1C=O(24352)',
    structure = SMILES('[CH2]C1([CH2])C=C1C=O'),
    E0 = (486.979,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
    molecularWeight = (94.1112,'amu'),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[1.46488,0.0571637,-4.50927e-05,1.74348e-08,-2.76561e-12,58659.6,22.6939], Tmin=(100,'K'), Tmax=(1438.65,'K')), NASAPolynomial(coeffs=[12.5302,0.0263977,-1.30146e-05,2.5698e-09,-1.82451e-13,55475.8,-34.7143], Tmin=(1438.65,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(486.979,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(295.164,'J/(mol*K)'), comment="""Thermo group additivity estimation: group(Cs-(Cds-Cds)(Cds-Cds)CsCs) + group(Cs-CsHHH) + group(Cs-CsHHH) + group(Cd-CdCs(CO)) + group(Cds-CdsCsH) + group(Cds-O2d(Cds-Cds)H) + ring(Cyclopropene) + radical(Neopentyl) + radical(Neopentyl)"""),
)

species(
    label = '[CH2]C(=C)C#CC=O(24384)',
    structure = SMILES('[CH2]C(=C)C#CC=O'),
    E0 = (251.79,'kJ/mol'),
    modes = [
        HarmonicOscillator(frequencies=([2950,3100,1380,975,1025,1650,2782.5,750,1395,475,1775,1000,350,440,435,1725,2100,2250,500,550,3000,3100,440,815,1455,1000,718.487],'cm^-1')),
        HinderedRotor(inertia=(0.336442,'amu*angstrom^2'), symmetry=1, barrier=(26.218,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(1.6034,'amu*angstrom^2'), symmetry=1, barrier=(36.8653,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(1.14011,'amu*angstrom^2'), symmetry=1, barrier=(26.2134,'kJ/mol'), semiclassical=False),
    ],
    spinMultiplicity = 2,
    opticalIsomers = 1,
    molecularWeight = (93.1033,'amu'),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[1.28996,0.0523709,-4.06463e-05,1.58218e-08,-2.46502e-12,30386.8,22.6167], Tmin=(100,'K'), Tmax=(1521.38,'K')), NASAPolynomial(coeffs=[13.9916,0.0189759,-7.72048e-06,1.39377e-09,-9.41302e-14,26522.1,-43.991], Tmin=(1521.38,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(251.79,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(270.22,'J/(mol*K)'), comment="""Thermo group additivity estimation: group(Cs-(Cds-Cds)HHH) + group(Cds-CdsCtCs) + group(Cds-CdsHH) + group(Cds-O2d(Cds-Cds)H) + group(Ct-Ct(Cds-Cds)) + group(Ct-CtCs) + radical(Allyl_P)"""),
)

species(
    label = '[CH2]C(=C)[C]=CC=O(24385)',
    structure = SMILES('[CH2]C(=C)[C]=CC=O'),
    E0 = (283.337,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
    molecularWeight = (94.1112,'amu'),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[1.00844,0.0650933,-6.12243e-05,2.99903e-08,-5.97224e-12,34185.9,22.2899], Tmin=(100,'K'), Tmax=(1195.27,'K')), NASAPolynomial(coeffs=[12.888,0.0253381,-1.13335e-05,2.16342e-09,-1.52033e-13,31346,-37.141], Tmin=(1195.27,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(283.337,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(295.164,'J/(mol*K)'), comment="""Thermo group additivity estimation: group(Cs-(Cds-Cds)HHH) + group(Cds-Cds(Cds-Cds)Cs) + group(Cds-Cds(Cds-Cds)H) + group(Cd-Cd(CO)H) + group(Cds-CdsHH) + group(Cds-O2d(Cds-Cds)H) + radical(Allyl_P) + radical(C=CJC=C)"""),
)

species(
    label = '[CH2]C([CH2])=CC=C=O(24386)',
    structure = SMILES('[CH2]C([CH2])=CC=C=O'),
    E0 = (189.445,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
    molecularWeight = (94.1112,'amu'),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[1.49714,0.0499096,-3.25917e-05,1.12628e-08,-1.6237e-12,22879.2,24.9779], Tmin=(100,'K'), Tmax=(1576.18,'K')), NASAPolynomial(coeffs=[10.8098,0.026276,-1.01004e-05,1.74983e-09,-1.1483e-13,19943.5,-24.1878], Tmin=(1576.18,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(189.445,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(295.164,'J/(mol*K)'), comment="""Thermo group additivity estimation: group(Cs-(Cds-Cds)HHH) + group(Cs-(Cds-Cds)HHH) + group(Cds-CdsCsCs) + group(Cd-Cd(CCO)H) + group(Cds-(Cdd-O2d)CsH) + radical(C=CC=CCJ) + radical(C=CC=CCJ)"""),
)

species(
    label = '[CH]C(=C)C=CC=O(20222)',
    structure = SMILES('[CH]C(=C)C=CC=O'),
    E0 = (303.527,'kJ/mol'),
    modes = [
        HarmonicOscillator(frequencies=([2950,3100,1380,975,1025,1650,2782.5,750,1395,475,1775,1000,350,440,435,1725,2995,3025,975,1000,1300,1375,400,500,1630,1680,180,180,180,180],'cm^-1')),
        HinderedRotor(inertia=(2.05674,'amu*angstrom^2'), symmetry=1, barrier=(47.2885,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(2.04662,'amu*angstrom^2'), symmetry=1, barrier=(47.0559,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(2.05042,'amu*angstrom^2'), symmetry=1, barrier=(47.1433,'kJ/mol'), semiclassical=False),
    ],
    spinMultiplicity = 3,
    opticalIsomers = 1,
    molecularWeight = (94.1112,'amu'),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[0.908717,0.0629631,-4.66694e-05,1.71276e-08,-2.54983e-12,36621,23.8842], Tmin=(100,'K'), Tmax=(1554.92,'K')), NASAPolynomial(coeffs=[14.9303,0.0268927,-1.18728e-05,2.20862e-09,-1.51161e-13,32260.5,-49.9514], Tmin=(1554.92,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(303.527,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(295.164,'J/(mol*K)'), comment="""Thermo group additivity estimation: group(Cs-(Cds-Cds)HHH) + group(Cds-Cds(Cds-Cds)Cs) + group(Cds-Cds(Cds-Cds)H) + group(Cd-Cd(CO)H) + group(Cds-CdsHH) + group(Cds-O2d(Cds-Cds)H) + radical(AllylJ2_triplet)"""),
)

species(
    label = 'C=C1C=C(C=O)C1(24297)',
    structure = SMILES('C=C1C=C(C=O)C1'),
    E0 = (88.9184,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
    molecularWeight = (94.1112,'amu'),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[1.55211,0.038041,2.02027e-05,-5.50418e-08,2.35949e-11,10796.5,19.2327], Tmin=(100,'K'), Tmax=(997.037,'K')), NASAPolynomial(coeffs=[15.7798,0.0187322,-7.57284e-06,1.52597e-09,-1.16196e-13,6082.07,-58.7801], Tmin=(997.037,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(88.9184,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(303.478,'J/(mol*K)'), comment="""Thermo group additivity estimation: group(Cs-(Cds-Cds)(Cds-Cds)HH) + group(Cd-CdCs(CO)) + group(Cds-Cds(Cds-Cds)Cs) + group(Cds-Cds(Cds-Cds)H) + group(Cds-O2d(Cds-Cds)H) + group(Cds-CdsHH) + ring(Cyclobutene)"""),
)

species(
    label = '[CH]=O(373)',
    structure = SMILES('[CH]=O'),
    E0 = (33.3613,'kJ/mol'),
    modes = [
        HarmonicOscillator(frequencies=([880.033,2045.55,4000],'cm^-1')),
    ],
    spinMultiplicity = 2,
    opticalIsomers = 1,
    molecularWeight = (29.018,'amu'),
    collisionModel = TransportData(shapeIndex=2, epsilon=(4140.62,'J/mol'), sigma=(3.59,'angstroms'), dipoleMoment=(0,'C*m'), polarizability=(0,'angstroms^3'), rotrelaxcollnum=0.0, comment="""GRI-Mech"""),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[4.35602,-0.00347087,1.25664e-05,-9.99488e-09,2.27888e-12,3995.77,2.75113], Tmin=(100,'K'), Tmax=(1565.72,'K')), NASAPolynomial(coeffs=[4.61867,0.00504456,-4.39241e-06,9.73282e-10,-7.07436e-14,2787.5,-2.22961], Tmin=(1565.72,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(33.3613,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(58.2013,'J/(mol*K)'), comment="""Thermo group additivity estimation: group(Cds-OdHH) + radical(HCdsJO)"""),
)

species(
    label = '[C]=CC([CH2])=C(17777)',
    structure = SMILES('[C]=CC([CH2])=C'),
    E0 = (765.932,'kJ/mol'),
    modes = [
        HarmonicOscillator(frequencies=([2950,3100,1380,975,1025,1650,3010,987.5,1337.5,450,1655,350,440,435,1725,3000,3100,440,815,1455,1000,180],'cm^-1')),
        HinderedRotor(inertia=(1.18271,'amu*angstrom^2'), symmetry=1, barrier=(27.1927,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(1.18148,'amu*angstrom^2'), symmetry=1, barrier=(27.1646,'kJ/mol'), semiclassical=False),
    ],
    spinMultiplicity = 4,
    opticalIsomers = 1,
    molecularWeight = (65.0932,'amu'),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[1.80667,0.0407524,-2.15729e-05,-5.89535e-09,6.53128e-12,92206,16.4335], Tmin=(100,'K'), Tmax=(963.658,'K')), NASAPolynomial(coeffs=[12.9302,0.0120127,-3.97215e-06,7.0025e-10,-4.97834e-14,89252.8,-41.019], Tmin=(963.658,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(765.932,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(224.491,'J/(mol*K)'), comment="""Thermo group additivity estimation: group(Cs-(Cds-Cds)HHH) + group(Cds-Cds(Cds-Cds)Cs) + group(Cds-Cds(Cds-Cds)H) + group(Cds-CdsHH) + group(Cds-CdsHH) + radical(Allyl_P) + radical(CdCdJ2_triplet)"""),
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
    E0 = (281.078,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS2',
    E0 = (429.158,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS3',
    E0 = (627.417,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS4',
    E0 = (539.901,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS5',
    E0 = (507.413,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS6',
    E0 = (474.68,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS7',
    E0 = (478.797,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS8',
    E0 = (505.99,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS9',
    E0 = (636.759,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS10',
    E0 = (342.271,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS11',
    E0 = (510.322,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS12',
    E0 = (665.339,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS13',
    E0 = (695.699,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS14',
    E0 = (712.526,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS15',
    E0 = (639.091,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS16',
    E0 = (385.256,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS17',
    E0 = (468.277,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS18',
    E0 = (411.523,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS19',
    E0 = (413.034,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS20',
    E0 = (412.017,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS21',
    E0 = (427.518,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS22',
    E0 = (306.052,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS23',
    E0 = (695.843,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS24',
    E0 = (409.479,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS25',
    E0 = (455.219,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS26',
    E0 = (669.138,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS27',
    E0 = (465.224,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS28',
    E0 = (626.435,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS29',
    E0 = (726.374,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS30',
    E0 = (533.279,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS31',
    E0 = (569.878,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS32',
    E0 = (781.851,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS33',
    E0 = (845.341,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS34',
    E0 = (486.979,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS35',
    E0 = (477.098,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS36',
    E0 = (496.234,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS37',
    E0 = (468.543,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS38',
    E0 = (443.858,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS39',
    E0 = (706.96,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS40',
    E0 = (289.363,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS41',
    E0 = (833.609,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

reaction(
    label = 'reaction1',
    reactants = ['[CH2]C(=C)C=C=C[O](22603)'],
    products = ['C=C=C(6077)', 'C#CC=O(21959)'],
    transitionState = 'TS1',
    kinetics = Arrhenius(A=(5e+12,'s^-1'), n=0, Ea=(0,'kJ/mol'), T0=(1,'K'), Tmin=(300,'K'), Tmax=(1500,'K'), comment="""Exact match found for rate rule [RJJ]
Euclidian distance = 0
family: 1,4_Linear_birad_scission"""),
)

reaction(
    label = 'reaction2',
    reactants = ['[CH2]C(=C)C=C=C[O](22603)'],
    products = ['C=C1CC1[C]=C[O](24299)'],
    transitionState = 'TS2',
    kinetics = Arrhenius(A=(1.264e+11,'s^-1'), n=0.35, Ea=(148.079,'kJ/mol'), T0=(1,'K'), comment="""Estimated using template [R4_S_D;doublebond_intra;radadd_intra_cs2H] for rate rule [R4_S_(Cd)_D;doublebond_intra;radadd_intra_cs2H]
Euclidian distance = 2.0
Multiplied by reaction path degeneracy 2.0
family: Intra_R_Add_Exocyclic
Ea raised from 147.0 to 148.1 kJ/mol to match endothermicity of reaction."""),
)

reaction(
    label = 'reaction3',
    reactants = ['[CH2]C(=C)C=C=C[O](22603)'],
    products = ['[CH2]C1([CH2])C=C=CO1(24365)'],
    transitionState = 'TS3',
    kinetics = Arrhenius(A=(1e+11,'s^-1'), n=0.21, Ea=(346.338,'kJ/mol'), T0=(1,'K'), comment="""Estimated using template [R6_SMM;doublebond_intra_2H;radadd_intra] for rate rule [R6_SMM;doublebond_intra_2H;radadd_intra_O]
Euclidian distance = 1.0
family: Intra_R_Add_Exocyclic"""),
)

reaction(
    label = 'reaction4',
    reactants = ['H(8)', '[CH2]C(=C)C=C=C=O(24366)'],
    products = ['[CH2]C(=C)C=C=C[O](22603)'],
    transitionState = 'TS4',
    kinetics = Arrhenius(A=(7048.5,'m^3/(mol*s)'), n=1.43107, Ea=(2.92155,'kJ/mol'), T0=(1,'K'), Tmin=(303.03,'K'), Tmax=(2000,'K'), comment="""From training reaction 188 used for Ck_O;HJ
Exact match found for rate rule [Ck_O;HJ]
Euclidian distance = 0
family: R_Addition_MultipleBond"""),
)

reaction(
    label = 'reaction5',
    reactants = ['C=C=C(6077)', '[CH]=C=C[O](8556)'],
    products = ['[CH2]C(=C)C=C=C[O](22603)'],
    transitionState = 'TS5',
    kinetics = Arrhenius(A=(1.04,'cm^3/(mol*s)'), n=3.05, Ea=(54.8104,'kJ/mol'), T0=(1,'K'), Tmin=(300,'K'), Tmax=(1500,'K'), comment="""From training reaction 2838 used for Ca_Cds-HH;CdsJ=Cdd
Exact match found for rate rule [Ca_Cds-HH;CdsJ=Cdd]
Euclidian distance = 0
family: R_Addition_MultipleBond"""),
)

reaction(
    label = 'reaction6',
    reactants = ['[CH2]C([CH2])=CC#CO(24367)'],
    products = ['[CH2]C(=C)C=C=C[O](22603)'],
    transitionState = 'TS6',
    kinetics = Arrhenius(A=(3.9024e+09,'s^-1'), n=1.09833, Ea=(177.686,'kJ/mol'), T0=(1,'K'), comment="""Estimated using template [R2H_S;Cd_rad_out_double;XH_out] for rate rule [R2H_S;Cd_rad_out_double;O_H_out]
Euclidian distance = 1.0
family: intra_H_migration"""),
)

reaction(
    label = 'reaction7',
    reactants = ['[CH2][C](C)C#CC=O(24368)'],
    products = ['[CH2]C(=C)C=C=C[O](22603)'],
    transitionState = 'TS7',
    kinetics = Arrhenius(A=(2.304e+09,'s^-1'), n=1.24, Ea=(151.879,'kJ/mol'), T0=(1,'K'), comment="""Estimated using template [R3H_SS;Cd_rad_out_double;Cs_H_out_2H] for rate rule [R3H_SS_2Cd;Cd_rad_out_double;Cs_H_out_2H]
Euclidian distance = 1.0
Multiplied by reaction path degeneracy 3.0
family: intra_H_migration"""),
)

reaction(
    label = 'reaction8',
    reactants = ['[CH]=C(C)C=C=C[O](24369)'],
    products = ['[CH2]C(=C)C=C=C[O](22603)'],
    transitionState = 'TS8',
    kinetics = Arrhenius(A=(720,'s^-1'), n=2.932, Ea=(129.315,'kJ/mol'), T0=(1,'K'), comment="""From training reaction 57 used for R3H_DS;Cd_rad_out_singleH;Cs_H_out_2H
Exact match found for rate rule [R3H_DS;Cd_rad_out_singleH;Cs_H_out_2H]
Euclidian distance = 0
Multiplied by reaction path degeneracy 3.0
family: intra_H_migration"""),
)

reaction(
    label = 'reaction9',
    reactants = ['[CH2]C([CH2])=C=C=CO(24370)'],
    products = ['[CH2]C(=C)C=C=C[O](22603)'],
    transitionState = 'TS9',
    kinetics = Arrhenius(A=(8.63625e+10,'s^-1'), n=1.0925, Ea=(294.328,'kJ/mol'), T0=(1,'K'), comment="""Estimated using template [R4H;Cd_rad_out_singleDe_Cd;XH_out] for rate rule [R4H_MMS;Cd_rad_out_singleDe_Cd;O_H_out]
Euclidian distance = 1.41421356237
family: intra_H_migration"""),
)

reaction(
    label = 'reaction10',
    reactants = ['[CH2]C(C)=C[C]=C=O(24371)'],
    products = ['[CH2]C(=C)C=C=C[O](22603)'],
    transitionState = 'TS10',
    kinetics = Arrhenius(A=(408000,'s^-1'), n=1.9199, Ea=(33.0402,'kJ/mol'), T0=(1,'K'), comment="""Estimated using an average for rate rule [R5H;Cd_rad_out;Cs_H_out_2H]
Euclidian distance = 0
Multiplied by reaction path degeneracy 3.0
family: intra_H_migration"""),
)

reaction(
    label = 'reaction11',
    reactants = ['[CH]C(=C)C=C=CO(24372)'],
    products = ['[CH2]C(=C)C=C=C[O](22603)'],
    transitionState = 'TS11',
    kinetics = Arrhenius(A=(1.86943e+06,'s^-1'), n=1.85754, Ea=(151.521,'kJ/mol'), T0=(1,'K'), comment="""Estimated using template [RnH;Cd_rad_out_singleH;XH_out] for rate rule [R6H;Cd_rad_out_singleH;O_H_out]
Euclidian distance = 2.2360679775
family: intra_H_migration"""),
)

reaction(
    label = 'reaction12',
    reactants = ['[CH2][C]=C(6078)', '[CH]=C=C[O](8556)'],
    products = ['[CH2]C(=C)C=C=C[O](22603)'],
    transitionState = 'TS12',
    kinetics = Arrhenius(A=(40577.9,'m^3/(mol*s)'), n=0.702818, Ea=(0,'kJ/mol'), T0=(1,'K'), comment="""Estimated using an average for rate rule [Cd_allenic;Y_rad]
Euclidian distance = 0
family: R_Recombination
Ea raised from -11.4 to 0 kJ/mol."""),
)

reaction(
    label = 'reaction13',
    reactants = ['H(8)', '[CH2]C([CH2])=C=C=C[O](24373)'],
    products = ['[CH2]C(=C)C=C=C[O](22603)'],
    transitionState = 'TS13',
    kinetics = Arrhenius(A=(2.048e+13,'cm^3/(mol*s)'), n=0.206, Ea=(0,'kJ/mol'), T0=(1,'K'), comment="""From training reaction 44 used for Cd_allenic;H_rad
Exact match found for rate rule [Cd_allenic;H_rad]
Euclidian distance = 0
family: R_Recombination
Ea raised from -0.7 to 0 kJ/mol."""),
)

reaction(
    label = 'reaction14',
    reactants = ['H(8)', '[CH]C(=C)C=C=C[O](22618)'],
    products = ['[CH2]C(=C)C=C=C[O](22603)'],
    transitionState = 'TS14',
    kinetics = Arrhenius(A=(5.78711e+07,'m^3/(mol*s)'), n=0.0433333, Ea=(0.458029,'kJ/mol'), T0=(1,'K'), comment="""Estimated using an average for rate rule [Cd_rad;H_rad]
Euclidian distance = 0
family: R_Recombination"""),
)

reaction(
    label = 'reaction15',
    reactants = ['H(8)', '[CH2]C([CH2])=C[C]=C=O(24374)'],
    products = ['[CH2]C(=C)C=C=C[O](22603)'],
    transitionState = 'TS15',
    kinetics = Arrhenius(A=(4.34078e+06,'m^3/(mol*s)'), n=0.278577, Ea=(0,'kJ/mol'), T0=(1,'K'), comment="""Estimated using an average for rate rule [H_rad;Y_rad]
Euclidian distance = 0
family: R_Recombination
Ea raised from -1.4 to 0 kJ/mol."""),
)

reaction(
    label = 'reaction16',
    reactants = ['[CH2]C(=C)C=C=C[O](22603)'],
    products = ['[O]C=C=C[C]1CC1(24375)'],
    transitionState = 'TS16',
    kinetics = Arrhenius(A=(2.03537e+13,'s^-1'), n=-0.296394, Ea=(104.177,'kJ/mol'), T0=(1,'K'), Tmin=(303.03,'K'), Tmax=(2000,'K'), comment="""Estimated using template [R3_D;doublebond_intra;radadd_intra_cs2H] for rate rule [R3_D;doublebond_intra_secDe_2H;radadd_intra_cs2H]
Euclidian distance = 2.0
family: Intra_R_Add_Endocyclic"""),
)

reaction(
    label = 'reaction17',
    reactants = ['[CH2]C(=C)C=C=C[O](22603)'],
    products = ['[CH2]C(=C)C=C1[CH]O1(24376)'],
    transitionState = 'TS17',
    kinetics = Arrhenius(A=(6.59543e+11,'s^-1'), n=0.400725, Ea=(187.199,'kJ/mol'), T0=(1,'K'), comment="""Estimated using template [R3_D;doublebond_intra;radadd_intra] for rate rule [R3_D;doublebond_intra;radadd_intra_O]
Euclidian distance = 1.0
family: Intra_R_Add_Endocyclic"""),
)

reaction(
    label = 'reaction18',
    reactants = ['[CH2]C(=C)C=C=C[O](22603)'],
    products = ['C=C1[CH]C(=C[O])C1(24377)'],
    transitionState = 'TS18',
    kinetics = Arrhenius(A=(6.43734e+08,'s^-1'), n=0.926191, Ea=(130.445,'kJ/mol'), T0=(1,'K'), comment="""Estimated using an average for rate rule [R4_S_D;doublebond_intra;radadd_intra_cs2H]
Euclidian distance = 0
Multiplied by reaction path degeneracy 2.0
family: Intra_R_Add_Endocyclic"""),
)

reaction(
    label = 'reaction19',
    reactants = ['[CH2]C(=C)C=C=C[O](22603)'],
    products = ['[CH2]C(=C)C1[C]=CO1(24378)'],
    transitionState = 'TS19',
    kinetics = Arrhenius(A=(2.47276e+11,'s^-1'), n=0.183155, Ea=(131.956,'kJ/mol'), T0=(1,'K'), comment="""Estimated using average of templates [Rn;doublebond_intra_CdCdd;radadd_intra] + [R4;doublebond_intra;radadd_intra] for rate rule [R4;doublebond_intra_CdCdd;radadd_intra_O]
Euclidian distance = 1.41421356237
family: Intra_R_Add_Endocyclic
Ea raised from 129.7 to 132.0 kJ/mol to match endothermicity of reaction."""),
)

reaction(
    label = 'reaction20',
    reactants = ['[CH2]C(=C)C=C=C[O](22603)'],
    products = ['C=C1C=[C]C([O])C1(24348)'],
    transitionState = 'TS20',
    kinetics = Arrhenius(A=(8.98919e+10,'s^-1'), n=0.314866, Ea=(130.939,'kJ/mol'), T0=(1,'K'), comment="""Estimated using template [R5_linear;doublebond_intra;radadd_intra_cs2H] for rate rule [R5_linear;doublebond_intra_CdCdd;radadd_intra_cs2H]
Euclidian distance = 1.0
Multiplied by reaction path degeneracy 2.0
family: Intra_R_Add_Endocyclic
Ea raised from 127.0 to 130.9 kJ/mol to match endothermicity of reaction."""),
)

reaction(
    label = 'reaction21',
    reactants = ['[CH2]C(=C)C=C=C[O](22603)'],
    products = ['[CH2]C1=C[C]=COC1(24334)'],
    transitionState = 'TS21',
    kinetics = Arrhenius(A=(2e+10,'s^-1'), n=0.19, Ea=(146.44,'kJ/mol'), T0=(1,'K'), comment="""Estimated using template [R6_SMM;doublebond_intra;radadd_intra] for rate rule [R6_SMM;doublebond_intra;radadd_intra_O]
Euclidian distance = 1.0
Multiplied by reaction path degeneracy 2.0
family: Intra_R_Add_Endocyclic"""),
)

reaction(
    label = 'reaction22',
    reactants = ['[CH2]C(=C)C=C=C[O](22603)'],
    products = ['C=C(C)C=C=C=O(24379)'],
    transitionState = 'TS22',
    kinetics = Arrhenius(A=(4.25221e+09,'s^-1'), n=0.137, Ea=(24.9733,'kJ/mol'), T0=(1,'K'), comment="""Estimated using template [R5;Y_rad;XH_Rrad] for rate rule [R5radExo;Y_rad;XH_Rrad]
Euclidian distance = 1.0
Multiplied by reaction path degeneracy 2.0
family: Intra_Disproportionation"""),
)

reaction(
    label = 'reaction23',
    reactants = ['[CH2]C([CH2])[C]=C=C[O](23592)'],
    products = ['[CH2]C(=C)C=C=C[O](22603)'],
    transitionState = 'TS23',
    kinetics = Arrhenius(A=(1.949e+11,'s^-1'), n=0.486, Ea=(22.8614,'kJ/mol'), T0=(1,'K'), comment="""Estimated using an average for rate rule [R2radExo;Y_rad_De;XH_Rrad]
Euclidian distance = 0
family: Intra_Disproportionation"""),
)

reaction(
    label = 'reaction24',
    reactants = ['[CH2]C([CH2])=[C]C=C[O](24380)'],
    products = ['[CH2]C(=C)C=C=C[O](22603)'],
    transitionState = 'TS24',
    kinetics = Arrhenius(A=(1.949e+11,'s^-1'), n=0.486, Ea=(22.8614,'kJ/mol'), T0=(1,'K'), comment="""Estimated using an average for rate rule [R2radExo;Y_rad_De;XH_Rrad]
Euclidian distance = 0
family: Intra_Disproportionation"""),
)

reaction(
    label = 'reaction25',
    reactants = ['[CH2]C([CH2])=CC=[C][O](24381)'],
    products = ['[CH2]C(=C)C=C=C[O](22603)'],
    transitionState = 'TS25',
    kinetics = Arrhenius(A=(1.4733e+10,'s^-1'), n=0.2847, Ea=(27.8529,'kJ/mol'), T0=(1,'K'), comment="""Estimated using average of templates [Rn;Y_rad;XH_Rrad_De] + [R2radExo;Y_rad;XH_Rrad] for rate rule [R2radExo;Y_rad;XH_Rrad_De]
Euclidian distance = 1.0
family: Intra_Disproportionation"""),
)

reaction(
    label = 'reaction26',
    reactants = ['[CH]C(=C)C[C]=C[O](20209)'],
    products = ['[CH2]C(=C)C=C=C[O](22603)'],
    transitionState = 'TS26',
    kinetics = Arrhenius(A=(4.00798e+09,'s^-1'), n=0.37, Ea=(78.2471,'kJ/mol'), T0=(1,'K'), comment="""Estimated using average of templates [R3;Y_rad;XH_Rrad_De] + [R3radExo;Y_rad;XH_Rrad] for rate rule [R3radExo;Y_rad;XH_Rrad_De]
Euclidian distance = 1.0
Multiplied by reaction path degeneracy 2.0
family: Intra_Disproportionation"""),
)

reaction(
    label = 'reaction27',
    reactants = ['[CH]C([CH2])=CC=C[O](20208)'],
    products = ['[CH2]C(=C)C=C=C[O](22603)'],
    transitionState = 'TS27',
    kinetics = Arrhenius(A=(5.14222e+08,'s^-1'), n=0.311, Ea=(24.9733,'kJ/mol'), T0=(1,'K'), comment="""Estimated using template [R4;Y_rad;XH_Rrad] for rate rule [R4radEndo;Y_rad;XH_Rrad]
Euclidian distance = 1.0
family: Intra_Disproportionation"""),
)

reaction(
    label = 'reaction28',
    reactants = ['[CH2]C([CH2])[CH][C]=C=O(23593)'],
    products = ['[CH2]C(=C)C=C=C[O](22603)'],
    transitionState = 'TS28',
    kinetics = Arrhenius(A=(5.14222e+08,'s^-1'), n=0.311, Ea=(24.9733,'kJ/mol'), T0=(1,'K'), comment="""Estimated using template [R4;Y_rad;XH_Rrad] for rate rule [R4radExo;Y_rad;XH_Rrad]
Euclidian distance = 1.0
family: Intra_Disproportionation"""),
)

reaction(
    label = 'reaction29',
    reactants = ['[CH]C(=C)C=[C]C[O](24382)'],
    products = ['[CH2]C(=C)C=C=C[O](22603)'],
    transitionState = 'TS29',
    kinetics = Arrhenius(A=(4.25221e+09,'s^-1'), n=0.137, Ea=(24.9733,'kJ/mol'), T0=(1,'K'), comment="""Estimated using template [R5;Y_rad;XH_Rrad] for rate rule [R5radEndo;Y_rad;XH_Rrad]
Euclidian distance = 1.0
Multiplied by reaction path degeneracy 2.0
family: Intra_Disproportionation"""),
)

reaction(
    label = 'reaction30',
    reactants = ['[CH2]C(C)=C[C]=[C][O](24383)'],
    products = ['[CH2]C(=C)C=C=C[O](22603)'],
    transitionState = 'TS30',
    kinetics = Arrhenius(A=(6.37831e+09,'s^-1'), n=0.137, Ea=(24.9733,'kJ/mol'), T0=(1,'K'), comment="""Estimated using template [R5;Y_rad;XH_Rrad] for rate rule [R5radEndo;Y_rad;XH_Rrad]
Euclidian distance = 1.0
Multiplied by reaction path degeneracy 3.0
family: Intra_Disproportionation"""),
)

reaction(
    label = 'reaction31',
    reactants = ['C=[C]CC=C=C[O](22607)'],
    products = ['[CH2]C(=C)C=C=C[O](22603)'],
    transitionState = 'TS31',
    kinetics = Arrhenius(A=(1.74842e+09,'s^-1'), n=1.084, Ea=(170.038,'kJ/mol'), T0=(1,'K'), comment="""Estimated using average of templates [cCsCJ;CdsJ;C] + [cCs(-HH)CJ;CJ;C] for rate rule [cCs(-HH)CJ;CdsJ;C]
Euclidian distance = 1.0
family: 1,2_shiftC"""),
)

reaction(
    label = 'reaction32',
    reactants = ['CH2(T)(28)', 'C=[C]C=C=C[O](23378)'],
    products = ['[CH2]C(=C)C=C=C[O](22603)'],
    transitionState = 'TS32',
    kinetics = Arrhenius(A=(1.14854e+06,'m^3/(mol*s)'), n=0.575199, Ea=(34.3157,'kJ/mol'), T0=(1,'K'), comment="""Estimated using template [Y_rad;Birad] for rate rule [Cd_rad/OneDe;Birad]
Euclidian distance = 3.0
family: Birad_R_Recombination"""),
)

reaction(
    label = 'reaction33',
    reactants = ['O(T)(63)', 'C#CC=C([CH2])[CH2](19657)'],
    products = ['[CH2]C(=C)C=C=C[O](22603)'],
    transitionState = 'TS33',
    kinetics = Arrhenius(A=(93609.6,'m^3/(mol*s)'), n=1.13083, Ea=(163.847,'kJ/mol'), T0=(1,'K'), Tmin=(303.03,'K'), Tmax=(2000,'K'), comment="""Estimated using template [Y_rad;O_birad] for rate rule [Cd_pri_rad;O_birad]
Euclidian distance = 2.0
family: Birad_R_Recombination"""),
)

reaction(
    label = 'reaction34',
    reactants = ['[CH2]C(=C)C=C=C[O](22603)'],
    products = ['[CH2]C1([CH2])C=C1C=O(24352)'],
    transitionState = 'TS34',
    kinetics = Arrhenius(A=(1.126e+14,'s^-1'), n=-0.355, Ea=(205.901,'kJ/mol'), T0=(1,'K'), comment="""Estimated using an average for rate rule [R4_D_D;doublebond_intra_2H;radadd_intra_cdsingleDe]
Euclidian distance = 0
family: Intra_R_Add_Exocyclic
Ea raised from 205.8 to 205.9 kJ/mol to match endothermicity of reaction."""),
)

reaction(
    label = 'reaction35',
    reactants = ['H(8)', '[CH2]C(=C)C#CC=O(24384)'],
    products = ['[CH2]C(=C)C=C=C[O](22603)'],
    transitionState = 'TS35',
    kinetics = Arrhenius(A=(2371.94,'m^3/(mol*s)'), n=1.49517, Ea=(13.5032,'kJ/mol'), T0=(1,'K'), comment="""Estimated using an average for rate rule [Ct-De_Ct-De;HJ]
Euclidian distance = 0
family: R_Addition_MultipleBond"""),
)

reaction(
    label = 'reaction36',
    reactants = ['[CH2][C]=C(6078)', 'C#CC=O(21959)'],
    products = ['[CH2]C(=C)C=C=C[O](22603)'],
    transitionState = 'TS36',
    kinetics = Arrhenius(A=(0.106247,'m^3/(mol*s)'), n=2.32278, Ea=(16.475,'kJ/mol'), T0=(1,'K'), comment="""Estimated using template [Ct-H_Ct-De;CJ] for rate rule [Ct-H_Ct-CO;CJ]
Euclidian distance = 1.0
family: R_Addition_MultipleBond"""),
)

reaction(
    label = 'reaction37',
    reactants = ['[CH2]C(=C)C=C=C[O](22603)'],
    products = ['[CH2]C(=C)[C]=CC=O(24385)'],
    transitionState = 'TS37',
    kinetics = Arrhenius(A=(7.7652e+07,'s^-1'), n=1.65613, Ea=(187.464,'kJ/mol'), T0=(1,'K'), comment="""Estimated using average of templates [R2H_D;Cd_rad_out_single;Cd_H_out_singleDe] + [R2H_D;Cd_rad_out_singleDe;Cd_H_out_single] for rate rule [R2H_D;Cd_rad_out_singleDe;Cd_H_out_singleDe]
Euclidian distance = 1.0
family: intra_H_migration"""),
)

reaction(
    label = 'reaction38',
    reactants = ['[CH2]C(=C)C=C=C[O](22603)'],
    products = ['[CH2]C([CH2])=CC=C=O(24386)'],
    transitionState = 'TS38',
    kinetics = Arrhenius(A=(4.96519e+09,'s^-1'), n=1.05826, Ea=(162.779,'kJ/mol'), T0=(1,'K'), comment="""Estimated using template [R2H_S;Cd_rad_out_Cd;XH_out] for rate rule [R2H_S;Cd_rad_out_Cd;CO_H_out]
Euclidian distance = 1.0
family: intra_H_migration"""),
)

reaction(
    label = 'reaction39',
    reactants = ['[CH]C(=C)C=CC=O(20222)'],
    products = ['[CH2]C(=C)C=C=C[O](22603)'],
    transitionState = 'TS39',
    kinetics = Arrhenius(A=(1.13764e+12,'s^-1'), n=1.09983, Ea=(403.433,'kJ/mol'), T0=(1,'K'), comment="""Estimated using template [R4H_DSD;Cd_rad_out_single;Cd_H_out_single] for rate rule [R4H_DSD;Cd_rad_out_singleH;Cd_H_out_singleDe]
Euclidian distance = 1.41421356237
family: intra_H_migration"""),
)

reaction(
    label = 'reaction40',
    reactants = ['[CH2]C(=C)C=C=C[O](22603)'],
    products = ['C=C1C=C(C=O)C1(24297)'],
    transitionState = 'TS40',
    kinetics = Arrhenius(A=(3.24e+12,'s^-1'), n=-0.305, Ea=(8.28432,'kJ/mol'), T0=(1,'K'), comment="""Estimated using template [R4;C_rad_out_2H;Ypri_rad_out] for rate rule [R4_SSD;C_rad_out_2H;CdsinglepriDe_rad_out]
Euclidian distance = 2.2360679775
Multiplied by reaction path degeneracy 2.0
family: Birad_recombination"""),
)

reaction(
    label = 'reaction41',
    reactants = ['[CH]=O(373)', '[C]=CC([CH2])=C(17777)'],
    products = ['[CH2]C(=C)C=C=C[O](22603)'],
    transitionState = 'TS41',
    kinetics = Arrhenius(A=(1.14854e+06,'m^3/(mol*s)'), n=0.575199, Ea=(34.3157,'kJ/mol'), T0=(1,'K'), comment="""Estimated using template [Y_rad;Birad] for rate rule [CO_pri_rad;Birad]
Euclidian distance = 2.0
family: Birad_R_Recombination"""),
)

network(
    label = '4905',
    isomers = [
        '[CH2]C(=C)C=C=C[O](22603)',
    ],
    reactants = [
        ('C=C=C(6077)', 'C#CC=O(21959)'),
    ],
    bathGas = {
        'N2': 0.25,
        'Ne': 0.25,
        'He': 0.25,
        'Ar': 0.25,
    },
)

pressureDependence(
    label = '4905',
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

