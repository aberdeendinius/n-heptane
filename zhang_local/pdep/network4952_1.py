species(
    label = '[CH]=C=COC=C=C[O](22649)',
    structure = SMILES('[CH]=C=COC=C=C[O]'),
    E0 = (338.118,'kJ/mol'),
    modes = [
        HarmonicOscillator(frequencies=([3120,650,792.5,1650,2995,3010,3025,975,987.5,1000,1300,1337.5,1375,400,450,500,1630,1655,1680,540,563.333,586.667,610,1970,2140,180,180,180],'cm^-1')),
        HinderedRotor(inertia=(1.5042,'amu*angstrom^2'), symmetry=1, barrier=(34.5844,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(1.52359,'amu*angstrom^2'), symmetry=1, barrier=(35.0304,'kJ/mol'), semiclassical=False),
    ],
    spinMultiplicity = 3,
    opticalIsomers = 1,
    molecularWeight = (108.095,'amu'),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[-1.00261,0.0837508,-9.25215e-05,4.69213e-08,-8.65362e-12,40868.9,32.6211], Tmin=(100,'K'), Tmax=(1585.25,'K')), NASAPolynomial(coeffs=[24.0008,0.00172769,3.00543e-06,-7.85906e-10,5.8134e-14,35320.5,-92.0229], Tmin=(1585.25,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(338.118,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(274.378,'J/(mol*K)'), comment="""Thermo group additivity estimation: group(O2s-(Cds-Cd)(Cds-Cd)) + group(O2s-(Cds-Cd)H) + group(Cds-CdsOsH) + group(Cds-CdsOsH) + group(Cds-CdsOsH) + group(Cds-CdsHH) + group(Cdd-CdsCds) + group(Cdd-CdsCds) + radical(C=COJ) + radical(C=C=CJ)"""),
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
    label = '[O]C=[C]C1C=C=CO1(25116)',
    structure = SMILES('[O]C=[C]C1C=C=CO1'),
    E0 = (524.68,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
    molecularWeight = (108.095,'amu'),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[0.731237,0.0470005,2.80629e-05,-9.16092e-08,4.41759e-11,63245.1,23.1692], Tmin=(100,'K'), Tmax=(931.324,'K')), NASAPolynomial(coeffs=[26.9262,-0.00172505,3.8147e-06,-7.17773e-10,3.83377e-14,55599.9,-116.193], Tmin=(931.324,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(524.68,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(278.535,'J/(mol*K)'), comment="""Thermo group additivity estimation: group(O2s-Cs(Cds-Cd)) + group(O2s-(Cds-Cd)H) + group(Cs-(Cds-Cds)(Cds-Cds)OsH) + group(Cds-CdsCsH) + group(Cds-CdsCsH) + group(Cds-CdsOsH) + group(Cds-CdsOsH) + group(Cdd-CdsCds) + ring(1,2-Cyclopentadiene) + radical(Cds_S) + radical(C=COJ)"""),
)

species(
    label = '[CH]=[C]C1OC=C=CO1(25171)',
    structure = SMILES('[CH]=[C]C1OC=C=CO1'),
    E0 = (402.049,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
    molecularWeight = (108.095,'amu'),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[-2.11481,0.125549,-0.000161709,9.15957e-08,-1.84297e-11,48583.7,12.6312], Tmin=(100,'K'), Tmax=(949.961,'K')), NASAPolynomial(coeffs=[31.3572,0.00404658,-5.47936e-07,3.55558e-11,-2.56782e-15,41347.2,-151.75], Tmin=(949.961,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(402.049,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(278.535,'J/(mol*K)'), comment="""Thermo group additivity estimation: group(O2s-Cs(Cds-Cd)) + group(O2s-Cs(Cds-Cd)) + group(Cs-CsOsOsH) + group(Cds-CdsCsH) + group(Cds-CdsOsH) + group(Cds-CdsOsH) + group(Cds-CdsHH) + group(Cdd-CdsCds) + ring(six-inringtwodouble-12) + radical(Cds_P) + radical(Cds_S)"""),
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
    label = '[CH]=C=COC=C=C=O(25172)',
    structure = SMILES('[CH]=C=COC=C=C=O'),
    E0 = (382.214,'kJ/mol'),
    modes = [
        HarmonicOscillator(frequencies=([3120,650,792.5,1650,2120,512.5,787.5,2995,3025,975,1000,1300,1375,400,500,1630,1680,540,563.333,586.667,610,1970,2140,180,4000],'cm^-1')),
        HinderedRotor(inertia=(0.160541,'amu*angstrom^2'), symmetry=1, barrier=(3.69116,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.160891,'amu*angstrom^2'), symmetry=1, barrier=(3.6992,'kJ/mol'), semiclassical=False),
    ],
    spinMultiplicity = 2,
    opticalIsomers = 1,
    molecularWeight = (107.087,'amu'),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[0.640224,0.0558998,-5.91583e-05,2.93136e-08,-5.30431e-12,46106.2,13.9268], Tmin=(100,'K'), Tmax=(1613,'K')), NASAPolynomial(coeffs=[17.0581,0.0029721,1.4198e-06,-4.18339e-10,3.14658e-14,42398.7,-68.2039], Tmin=(1613,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(382.214,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(249.434,'J/(mol*K)'), comment="""Thermo group additivity estimation: group(O2s-(Cds-Cd)(Cds-Cd)) + group(Cds-CdsOsH) + group(Cds-CdsOsH) + group(Cds-CdsHH) + group(Cdd-CdsCds) + group(Cdd-CdsCds) + radical(C=C=CJ)"""),
)

species(
    label = '[CH]=C=COC=C=[C]O(25173)',
    structure = SMILES('[CH]=C=COC=C=[C]O'),
    E0 = (436.399,'kJ/mol'),
    modes = [
        HarmonicOscillator(frequencies=([3120,650,792.5,1650,1685,370,540,563.333,586.667,610,1970,2140,3615,1277.5,1000,2995,3025,975,1000,1300,1375,400,500,1630,1680,180,180],'cm^-1')),
        HinderedRotor(inertia=(1.28678,'amu*angstrom^2'), symmetry=1, barrier=(29.5856,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(1.28783,'amu*angstrom^2'), symmetry=1, barrier=(29.6097,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(1.28693,'amu*angstrom^2'), symmetry=1, barrier=(29.5891,'kJ/mol'), semiclassical=False),
    ],
    spinMultiplicity = 3,
    opticalIsomers = 1,
    molecularWeight = (108.095,'amu'),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[-0.740437,0.0872955,-0.000105921,5.91989e-08,-1.21176e-11,52672.1,33.379], Tmin=(100,'K'), Tmax=(1395.68,'K')), NASAPolynomial(coeffs=[23.7782,0.00147679,3.02309e-06,-8.2172e-10,6.33433e-14,47342.4,-87.6584], Tmin=(1395.68,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(436.399,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(270.22,'J/(mol*K)'), comment="""Thermo group additivity estimation: group(O2s-(Cds-Cd)(Cds-Cd)) + group(O2s-(Cds-Cd)H) + group(Cds-CdsOsH) + group(Cds-CdsOsH) + group(Cds-CdsOsH) + group(Cds-CdsHH) + group(Cdd-CdsCds) + group(Cdd-CdsCds) + radical(C=CJO) + radical(C=C=CJ)"""),
)

species(
    label = '[CH2]C#COC=C=C[O](25174)',
    structure = SMILES('[CH2]C#COC=C=C[O]'),
    E0 = (366.615,'kJ/mol'),
    modes = [
        HarmonicOscillator(frequencies=([540,610,2055,2995,3025,975,1000,1300,1375,400,500,1630,1680,2100,2250,500,550,3000,3100,440,815,1455,1000,180,180,180,180],'cm^-1')),
        HinderedRotor(inertia=(1.43316,'amu*angstrom^2'), symmetry=1, barrier=(32.9512,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(1.43528,'amu*angstrom^2'), symmetry=1, barrier=(32.9999,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(1.43318,'amu*angstrom^2'), symmetry=1, barrier=(32.9517,'kJ/mol'), semiclassical=False),
    ],
    spinMultiplicity = 3,
    opticalIsomers = 1,
    molecularWeight = (108.095,'amu'),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[0.389937,0.0673063,-4.90794e-05,1.18222e-09,8.20468e-12,44234.6,28.0596], Tmin=(100,'K'), Tmax=(950.623,'K')), NASAPolynomial(coeffs=[21.2677,0.00849266,-2.09135e-06,3.59462e-10,-2.85907e-14,38953.3,-78.5066], Tmin=(950.623,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(366.615,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(270.22,'J/(mol*K)'), comment="""Thermo group additivity estimation: group(O2s-(Cds-Cd)(Cds-Cd)) + group(O2s-(Cds-Cd)H) + group(Cs-CtHHH) + group(Cds-CdsOsH) + group(Cds-CdsOsH) + group(Cdd-CdsCds) + group(Ct-CtCs) + group(Ct-CtOs) + radical(C=COJ) + radical(Propargyl)"""),
)

species(
    label = '[CH]=C=CO[C]=C=CO(25175)',
    structure = SMILES('[CH]=C=CO[C]=C=CO'),
    E0 = (436.399,'kJ/mol'),
    modes = [
        HarmonicOscillator(frequencies=([3120,650,792.5,1650,1685,370,540,563.333,586.667,610,1970,2140,3615,1277.5,1000,2995,3025,975,1000,1300,1375,400,500,1630,1680,180,180],'cm^-1')),
        HinderedRotor(inertia=(1.28678,'amu*angstrom^2'), symmetry=1, barrier=(29.5856,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(1.28783,'amu*angstrom^2'), symmetry=1, barrier=(29.6097,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(1.28693,'amu*angstrom^2'), symmetry=1, barrier=(29.5891,'kJ/mol'), semiclassical=False),
    ],
    spinMultiplicity = 3,
    opticalIsomers = 1,
    molecularWeight = (108.095,'amu'),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[-0.740437,0.0872955,-0.000105921,5.91989e-08,-1.21176e-11,52672.1,33.379], Tmin=(100,'K'), Tmax=(1395.68,'K')), NASAPolynomial(coeffs=[23.7782,0.00147679,3.02309e-06,-8.2172e-10,6.33433e-14,47342.4,-87.6584], Tmin=(1395.68,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(436.399,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(270.22,'J/(mol*K)'), comment="""Thermo group additivity estimation: group(O2s-(Cds-Cd)(Cds-Cd)) + group(O2s-(Cds-Cd)H) + group(Cds-CdsOsH) + group(Cds-CdsOsH) + group(Cds-CdsOsH) + group(Cds-CdsHH) + group(Cdd-CdsCds) + group(Cdd-CdsCds) + radical(C=CJO) + radical(C=C=CJ)"""),
)

species(
    label = 'C=[C][CH]OC#CC=O(25176)',
    structure = SMILES('C=[C][CH]OC#CC=O'),
    E0 = (352.199,'kJ/mol'),
    modes = [
        HarmonicOscillator(frequencies=([1685,370,3025,407.5,1350,352.5,2782.5,750,1395,475,1775,1000,2100,2250,500,550,2950,3100,1380,975,1025,1650,370.387,370.388,370.389,370.39],'cm^-1')),
        HinderedRotor(inertia=(0.283691,'amu*angstrom^2'), symmetry=1, barrier=(27.6176,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.283697,'amu*angstrom^2'), symmetry=1, barrier=(27.6176,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.283692,'amu*angstrom^2'), symmetry=1, barrier=(27.6176,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.28369,'amu*angstrom^2'), symmetry=1, barrier=(27.6176,'kJ/mol'), semiclassical=False),
    ],
    spinMultiplicity = 3,
    opticalIsomers = 1,
    molecularWeight = (108.095,'amu'),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[1.03346,0.0661528,-6.73599e-05,3.49807e-08,-7.33713e-12,42465.8,26.1318], Tmin=(100,'K'), Tmax=(1140.61,'K')), NASAPolynomial(coeffs=[13.1822,0.0235489,-1.13328e-05,2.23424e-09,-1.59836e-13,39694.4,-34.0776], Tmin=(1140.61,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(352.199,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(266.063,'J/(mol*K)'), comment="""Thermo group additivity estimation: group(O2s-Cs(Cds-Cd)) + group(Cs-(Cds-Cds)OsHH) + group(Cds-CdsCsH) + group(Cds-CdsHH) + group(Cds-O2d(Cds-Cds)H) + group(Ct-CtCs) + group(Ct-CtOs) + radical(C=CCJ(O)C) + radical(Cds_S)"""),
)

species(
    label = '[CH]=C=[C]OC=C=CO(25177)',
    structure = SMILES('[CH]=C=[C]OC=C=CO'),
    E0 = (436.399,'kJ/mol'),
    modes = [
        HarmonicOscillator(frequencies=([3120,650,792.5,1650,1685,370,540,563.333,586.667,610,1970,2140,3615,1277.5,1000,2995,3025,975,1000,1300,1375,400,500,1630,1680,180,180],'cm^-1')),
        HinderedRotor(inertia=(1.28245,'amu*angstrom^2'), symmetry=1, barrier=(29.4861,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(1.29091,'amu*angstrom^2'), symmetry=1, barrier=(29.6805,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(1.28376,'amu*angstrom^2'), symmetry=1, barrier=(29.5161,'kJ/mol'), semiclassical=False),
    ],
    spinMultiplicity = 3,
    opticalIsomers = 1,
    molecularWeight = (108.095,'amu'),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[-0.740446,0.0872956,-0.000105921,5.91992e-08,-1.21177e-11,52672.1,33.379], Tmin=(100,'K'), Tmax=(1395.67,'K')), NASAPolynomial(coeffs=[23.7784,0.00147654,3.02322e-06,-8.21748e-10,6.33456e-14,47342.3,-87.6593], Tmin=(1395.67,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(436.399,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(270.22,'J/(mol*K)'), comment="""Thermo group additivity estimation: group(O2s-(Cds-Cd)(Cds-Cd)) + group(O2s-(Cds-Cd)H) + group(Cds-CdsOsH) + group(Cds-CdsOsH) + group(Cds-CdsOsH) + group(Cds-CdsHH) + group(Cdd-CdsCds) + group(Cdd-CdsCds) + radical(C=C=CJ) + radical(C=CJO)"""),
)

species(
    label = 'C=C=CO[CH][C]=C=O(25178)',
    structure = SMILES('C=C=CO[CH][C]=C=O'),
    E0 = (322.897,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
    molecularWeight = (108.095,'amu'),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[0.226284,0.0822748,-0.000104567,6.61249e-08,-1.6241e-11,38972.2,26.6304], Tmin=(100,'K'), Tmax=(1004.05,'K')), NASAPolynomial(coeffs=[16.6478,0.0168549,-6.83499e-06,1.23413e-09,-8.40591e-14,35674.5,-52.6606], Tmin=(1004.05,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(322.897,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(270.22,'J/(mol*K)'), comment="""Thermo group additivity estimation: group(O2s-Cs(Cds-Cd)) + group(Cs-(Cds-Cdd-O2d)OsHH) + group(Cds-(Cdd-O2d)CsH) + group(Cds-CdsOsH) + group(Cds-CdsHH) + group(Cdd-CdsCds) + radical(C=CCJ(O)C) + radical(Cds_S)"""),
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
    label = '[O]C=C=C[O](22349)',
    structure = SMILES('[O]C=C=C[O]'),
    E0 = (48.0679,'kJ/mol'),
    modes = [
        HarmonicOscillator(frequencies=([2995,3025,975,1000,1300,1375,400,500,1630,1680,540,610,2055,180,180],'cm^-1')),
    ],
    spinMultiplicity = 3,
    opticalIsomers = 1,
    molecularWeight = (70.0468,'amu'),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[2.10387,0.0254095,2.29542e-05,-6.61705e-08,3.24893e-11,5864.8,15.5654], Tmin=(100,'K'), Tmax=(907.068,'K')), NASAPolynomial(coeffs=[19.4725,-0.00784929,6.29342e-06,-1.25746e-09,8.23938e-14,931.198,-76.3609], Tmin=(907.068,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(48.0679,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(157.975,'J/(mol*K)'), comment="""Thermo group additivity estimation: group(O2s-(Cds-Cd)H) + group(O2s-(Cds-Cd)H) + group(Cds-CdsOsH) + group(Cds-CdsOsH) + group(Cdd-CdsCds) + radical(C=COJ) + radical(C=COJ)"""),
)

species(
    label = '[CH]=C=CO[C]=C=C[O](25179)',
    structure = SMILES('[CH]=C=CO[C]=C=C[O]'),
    E0 = (577.862,'kJ/mol'),
    modes = [
        HarmonicOscillator(frequencies=([3120,650,792.5,1650,1685,370,2995,3025,975,1000,1300,1375,400,500,1630,1680,540,563.333,586.667,610,1970,2140,180,180,180],'cm^-1')),
        HinderedRotor(inertia=(1.47578,'amu*angstrom^2'), symmetry=1, barrier=(33.9311,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(1.48058,'amu*angstrom^2'), symmetry=1, barrier=(34.0415,'kJ/mol'), semiclassical=False),
    ],
    spinMultiplicity = 4,
    opticalIsomers = 1,
    molecularWeight = (107.087,'amu'),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[-0.0888606,0.077437,-9.21207e-05,5.14584e-08,-1.06847e-11,69658.7,32.4543], Tmin=(100,'K'), Tmax=(1336.83,'K')), NASAPolynomial(coeffs=[20.9797,0.00456396,6.79884e-07,-3.22491e-10,2.77176e-14,64904.3,-72.0194], Tmin=(1336.83,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(577.862,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(249.434,'J/(mol*K)'), comment="""Thermo group additivity estimation: group(O2s-(Cds-Cd)(Cds-Cd)) + group(O2s-(Cds-Cd)H) + group(Cds-CdsOsH) + group(Cds-CdsOsH) + group(Cds-CdsOsH) + group(Cds-CdsHH) + group(Cdd-CdsCds) + group(Cdd-CdsCds) + radical(C=C=CJ) + radical(C=CJO) + radical(C=COJ)"""),
)

species(
    label = '[CH]=C=[C]OC=C=C[O](25180)',
    structure = SMILES('[CH]=C=[C]OC=C=C[O]'),
    E0 = (577.862,'kJ/mol'),
    modes = [
        HarmonicOscillator(frequencies=([3120,650,792.5,1650,1685,370,2995,3025,975,1000,1300,1375,400,500,1630,1680,540,563.333,586.667,610,1970,2140,180,180,180],'cm^-1')),
        HinderedRotor(inertia=(1.4798,'amu*angstrom^2'), symmetry=1, barrier=(34.0235,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(1.4755,'amu*angstrom^2'), symmetry=1, barrier=(33.9246,'kJ/mol'), semiclassical=False),
    ],
    spinMultiplicity = 4,
    opticalIsomers = 1,
    molecularWeight = (107.087,'amu'),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[-0.088843,0.0774368,-9.21202e-05,5.14579e-08,-1.06845e-11,69658.7,32.4543], Tmin=(100,'K'), Tmax=(1336.85,'K')), NASAPolynomial(coeffs=[20.9794,0.00456435,6.79677e-07,-3.22445e-10,2.77139e-14,64904.5,-72.0179], Tmin=(1336.85,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(577.862,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(249.434,'J/(mol*K)'), comment="""Thermo group additivity estimation: group(O2s-(Cds-Cd)(Cds-Cd)) + group(O2s-(Cds-Cd)H) + group(Cds-CdsOsH) + group(Cds-CdsOsH) + group(Cds-CdsOsH) + group(Cds-CdsHH) + group(Cdd-CdsCds) + group(Cdd-CdsCds) + radical(C=COJ) + radical(C=CJO) + radical(C=C=CJ)"""),
)

species(
    label = '[CH]=C=CO[CH][C]=C=O(25181)',
    structure = SMILES('[CH]=C=CO[CH][C]=C=O'),
    E0 = (477.374,'kJ/mol'),
    modes = [
        HarmonicOscillator(frequencies=([3120,650,792.5,1650,3025,407.5,1350,352.5,540,610,2055,2120,512.5,787.5,1685,370,3010,987.5,1337.5,450,1655,180,180,180],'cm^-1')),
        HinderedRotor(inertia=(1.49853,'amu*angstrom^2'), symmetry=1, barrier=(34.4542,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(1.49137,'amu*angstrom^2'), symmetry=1, barrier=(34.2896,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(1.49516,'amu*angstrom^2'), symmetry=1, barrier=(34.3767,'kJ/mol'), semiclassical=False),
    ],
    spinMultiplicity = 4,
    opticalIsomers = 1,
    molecularWeight = (107.087,'amu'),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[0.0661567,0.0853807,-0.000120168,8.23023e-08,-2.15244e-11,57557.5,27.8267], Tmin=(100,'K'), Tmax=(980.803,'K')), NASAPolynomial(coeffs=[17.562,0.0115844,-3.57069e-06,5.09329e-10,-2.85159e-14,54243,-55.6424], Tmin=(980.803,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(477.374,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(245.277,'J/(mol*K)'), comment="""Thermo group additivity estimation: group(O2s-Cs(Cds-Cd)) + group(Cs-(Cds-Cdd-O2d)OsHH) + group(Cds-(Cdd-O2d)CsH) + group(Cds-CdsOsH) + group(Cds-CdsHH) + group(Cdd-CdsCds) + radical(C=C=CJ) + radical(C=CCJ(O)C) + radical(Cds_S)"""),
)

species(
    label = '[CH]=C=COC=C1[CH]O1(25182)',
    structure = SMILES('[CH]=C=COC=C1[CH]O1'),
    E0 = (364.3,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
    molecularWeight = (108.095,'amu'),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[0.543351,0.0494125,3.03201e-05,-1.04549e-07,5.25552e-11,43964.6,25.6453], Tmin=(100,'K'), Tmax=(904.532,'K')), NASAPolynomial(coeffs=[29.8575,-0.00900948,9.1131e-06,-1.88341e-09,1.24674e-13,35748.4,-128.94], Tmin=(904.532,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(364.3,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(274.378,'J/(mol*K)'), comment="""Thermo group additivity estimation: group(O2s-Cs(Cds-Cd)) + group(O2s-(Cds-Cd)(Cds-Cd)) + group(Cs-(Cds-Cds)OsHH) + group(Cds-CdsCsOs) + group(Cds-CdsOsH) + group(Cds-CdsOsH) + group(Cds-CdsHH) + group(Cdd-CdsCds) + ring(methyleneoxirane) + radical(C=C=CJ) + radical(C=CCJO)"""),
)

species(
    label = '[O]C=C=COC1[C]=C1(25183)',
    structure = SMILES('[O]C=C=COC1[C]=C1'),
    E0 = (469.349,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
    molecularWeight = (108.095,'amu'),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[0.152288,0.0637268,-1.63902e-05,-4.91938e-08,3.02101e-11,56607.6,26.4181], Tmin=(100,'K'), Tmax=(925.262,'K')), NASAPolynomial(coeffs=[28.1754,-0.00324373,4.35165e-06,-8.56918e-10,5.1433e-14,49102.8,-119.132], Tmin=(925.262,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(469.349,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(274.378,'J/(mol*K)'), comment="""Thermo group additivity estimation: group(O2s-Cs(Cds-Cd)) + group(O2s-(Cds-Cd)H) + group(Cs-(Cds-Cds)(Cds-Cds)OsH) + group(Cds-CdsCsH) + group(Cds-CdsCsH) + group(Cds-CdsOsH) + group(Cds-CdsOsH) + group(Cdd-CdsCds) + ring(Cyclopropene) + radical(cyclopropenyl-vinyl) + radical(C=COJ)"""),
)

species(
    label = '[CH]=C=COC1[C]=CO1(25184)',
    structure = SMILES('[CH]=C=COC1[C]=CO1'),
    E0 = (401.858,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
    molecularWeight = (108.095,'amu'),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[0.151652,0.0587124,9.90722e-06,-8.65939e-08,4.66445e-11,48495.2,24.1003], Tmin=(100,'K'), Tmax=(907.689,'K')), NASAPolynomial(coeffs=[31.5194,-0.0101692,9.13367e-06,-1.85322e-09,1.21686e-13,39943.9,-139.93], Tmin=(907.689,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(401.858,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(274.378,'J/(mol*K)'), comment="""Thermo group additivity estimation: group(O2s-Cs(Cds-Cd)) + group(O2s-Cs(Cds-Cd)) + group(Cs-CsOsOsH) + group(Cds-CdsCsH) + group(Cds-CdsOsH) + group(Cds-CdsOsH) + group(Cds-CdsHH) + group(Cdd-CdsCds) + ring(Cyclobutene) + radical(Cds_S) + radical(C=C=CJ)"""),
)

species(
    label = '[O]C=C1[CH]OC=C=C1(25114)',
    structure = SMILES('[O]C=C1[CH]OC=C=C1'),
    E0 = (117.112,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
    molecularWeight = (108.095,'amu'),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[1.22044,0.0273263,9.77896e-05,-1.64995e-07,6.95734e-11,14216.5,17.0539], Tmin=(100,'K'), Tmax=(940.714,'K')), NASAPolynomial(coeffs=[27.6409,0.000980083,2.6763e-06,-4.13076e-10,9.6179e-15,5440.61,-129.019], Tmin=(940.714,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(117.112,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(282.692,'J/(mol*K)'), comment="""Thermo group additivity estimation: group(O2s-Cs(Cds-Cd)) + group(O2s-(Cds-Cd)H) + group(Cs-(Cds-Cds)OsHH) + group(Cds-Cds(Cds-Cds)Cs) + group(Cds-Cds(Cds-Cds)H) + group(Cds-CdsOsH) + group(Cds-CdsOsH) + group(Cdd-CdsCds) + ring(Cyclohexane) + radical(C=COJ) + radical(C=CCJ(O)C)"""),
)

species(
    label = '[O]C1[C]=COC=C=C1(25185)',
    structure = SMILES('[O]C1[C]=COC=C=C1'),
    E0 = (487.062,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
    molecularWeight = (108.095,'amu'),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[1.38302,0.0422373,9.80236e-06,-4.81247e-08,2.24607e-11,58688,20.7611], Tmin=(100,'K'), Tmax=(975.204,'K')), NASAPolynomial(coeffs=[16.9168,0.015359,-5.51531e-06,1.0806e-09,-8.30382e-14,53906.6,-62.7712], Tmin=(975.204,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(487.062,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(282.692,'J/(mol*K)'), comment="""Thermo group additivity estimation: group(O2s-(Cds-Cd)(Cds-Cd)) + group(O2s-CsH) + group(Cs-(Cds-Cds)(Cds-Cds)OsH) + group(Cds-CdsCsH) + group(Cds-CdsCsH) + group(Cds-CdsOsH) + group(Cds-CdsOsH) + group(Cdd-CdsCds) + ring(1_2_cycloheptadiene) + radical(CC(C)OJ) + radical(Cds_S)"""),
)

species(
    label = '[CH]=C1[CH]OC=C=CO1(25157)',
    structure = SMILES('[CH]=C1[CH]OC=C=CO1'),
    E0 = (384.263,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
    molecularWeight = (108.095,'amu'),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[1.12454,0.0449135,1.47534e-05,-6.18771e-08,2.95353e-11,46336.4,19.5213], Tmin=(100,'K'), Tmax=(951.702,'K')), NASAPolynomial(coeffs=[20.1486,0.0104418,-2.60661e-06,5.0364e-10,-4.35242e-14,40655.5,-82.1392], Tmin=(951.702,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(384.263,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(282.692,'J/(mol*K)'), comment="""Thermo group additivity estimation: group(O2s-Cs(Cds-Cd)) + group(O2s-(Cds-Cd)(Cds-Cd)) + group(Cs-(Cds-Cds)OsHH) + group(Cds-CdsCsOs) + group(Cds-CdsOsH) + group(Cds-CdsOsH) + group(Cds-CdsHH) + group(Cdd-CdsCds) + ring(1_2_cycloheptadiene) + radical(Cds_P) + radical(C=CCJ(O)C)"""),
)

species(
    label = 'C=C=COC=C=C=O(25186)',
    structure = SMILES('C=C=COC=C=C=O'),
    E0 = (227.738,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
    molecularWeight = (108.095,'amu'),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[1.43135,0.0455878,-1.93431e-05,-1.70463e-08,1.23886e-11,27493,10.4497], Tmin=(100,'K'), Tmax=(948.713,'K')), NASAPolynomial(coeffs=[16.7996,0.00758923,-1.63345e-06,2.82389e-10,-2.35629e-14,23371,-69.24], Tmin=(948.713,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(227.738,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(274.378,'J/(mol*K)'), comment="""Thermo group additivity estimation: group(O2s-(Cds-Cd)(Cds-Cd)) + group(Cds-CdsOsH) + group(Cds-CdsOsH) + group(Cds-CdsHH) + group(Cdd-CdsCds) + group(Cdd-CdsCds)"""),
)

species(
    label = '[CH]=C=CO[C]C=C[O](25187)',
    structure = SMILES('[CH]=C=CO[C]C=C[O]'),
    E0 = (612.766,'kJ/mol'),
    modes = [
        HarmonicOscillator(frequencies=([3120,650,792.5,1650,2995,3010,3025,975,987.5,1000,1300,1337.5,1375,400,450,500,1630,1655,1680,540,610,2055,180,180,180,180,180],'cm^-1')),
        HinderedRotor(inertia=(1.31463,'amu*angstrom^2'), symmetry=1, barrier=(30.226,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(1.3132,'amu*angstrom^2'), symmetry=1, barrier=(30.1931,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(1.31249,'amu*angstrom^2'), symmetry=1, barrier=(30.1768,'kJ/mol'), semiclassical=False),
    ],
    spinMultiplicity = 5,
    opticalIsomers = 1,
    molecularWeight = (108.095,'amu'),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[-1.60809,0.0915485,-0.00010498,5.35388e-08,-9.81909e-12,73928,33.4185], Tmin=(100,'K'), Tmax=(1610.48,'K')), NASAPolynomial(coeffs=[27.2314,-0.00354788,5.44955e-06,-1.222e-09,8.6138e-14,67682.2,-110.01], Tmin=(1610.48,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(612.766,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(270.22,'J/(mol*K)'), comment="""Thermo group additivity estimation: group(O2s-Cs(Cds-Cd)) + group(O2s-(Cds-Cd)H) + group(Cs-(Cds-Cds)OsHH) + group(Cds-CdsCsH) + group(Cds-CdsOsH) + group(Cds-CdsOsH) + group(Cds-CdsHH) + group(Cdd-CdsCds) + radical(CH2_triplet) + radical(C=COJ) + radical(C=C=CJ)"""),
)

species(
    label = '[CH]=C=CO[CH]C=[C][O](25188)',
    structure = SMILES('[CH]=C=CO[CH]C=[C][O]'),
    E0 = (502.461,'kJ/mol'),
    modes = [
        HarmonicOscillator(frequencies=([1685,370,3025,407.5,1350,352.5,2995,3025,975,1000,1300,1375,400,500,1630,1680,3120,650,792.5,1650,540,610,2055,180,180,180,180],'cm^-1')),
        HinderedRotor(inertia=(1.5124,'amu*angstrom^2'), symmetry=1, barrier=(34.7732,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(1.5127,'amu*angstrom^2'), symmetry=1, barrier=(34.7799,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(1.51528,'amu*angstrom^2'), symmetry=1, barrier=(34.8392,'kJ/mol'), semiclassical=False),
    ],
    spinMultiplicity = 5,
    opticalIsomers = 1,
    molecularWeight = (108.095,'amu'),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[-0.22241,0.0763094,-8.19546e-05,4.19309e-08,-8.01402e-12,60598.3,33.7649], Tmin=(100,'K'), Tmax=(1458.23,'K')), NASAPolynomial(coeffs=[21.2076,0.00733882,-5.29663e-07,-8.51523e-11,1.05442e-14,55431.4,-73.9928], Tmin=(1458.23,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(502.461,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(270.22,'J/(mol*K)'), comment="""Thermo group additivity estimation: group(O2s-Cs(Cds-Cd)) + group(O2s-(Cds-Cd)H) + group(Cs-(Cds-Cds)OsHH) + group(Cds-CdsCsH) + group(Cds-CdsOsH) + group(Cds-CdsOsH) + group(Cds-CdsHH) + group(Cdd-CdsCds) + radical(C=COJ) + radical(C=CJO) + radical(C=CCJ(O)C) + radical(C=C=CJ)"""),
)

species(
    label = '[CH]=[C]CO[C]=C=C[O](25189)',
    structure = SMILES('[CH]=[C]CO[C]=C=C[O]'),
    E0 = (721.983,'kJ/mol'),
    modes = [
        HarmonicOscillator(frequencies=([3120,650,792.5,1650,2750,2850,1437.5,1250,1305,750,350,540,610,2055,3010,987.5,1337.5,450,1655,1670,1700,300,440,180,180,180,180],'cm^-1')),
        HinderedRotor(inertia=(1.03756,'amu*angstrom^2'), symmetry=1, barrier=(23.8556,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(1.03721,'amu*angstrom^2'), symmetry=1, barrier=(23.8476,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(1.03723,'amu*angstrom^2'), symmetry=1, barrier=(23.848,'kJ/mol'), semiclassical=False),
    ],
    spinMultiplicity = 5,
    opticalIsomers = 1,
    molecularWeight = (108.095,'amu'),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[0.101912,0.0781432,-8.40962e-05,3.82995e-08,-4.87671e-12,86981.7,31.9116], Tmin=(100,'K'), Tmax=(947.838,'K')), NASAPolynomial(coeffs=[20.6494,0.00897785,-2.40869e-06,3.7637e-10,-2.59685e-14,82298.4,-70.2754], Tmin=(947.838,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(721.983,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(270.22,'J/(mol*K)'), comment="""Thermo group additivity estimation: group(O2s-Cs(Cds-Cd)) + group(O2s-(Cds-Cd)H) + group(Cs-(Cds-Cds)OsHH) + group(Cds-CdsCsH) + group(Cds-CdsOsH) + group(Cds-CdsHH) + group(Cds-CdsOsH) + group(Cdd-CdsCds) + radical(C=COJ) + radical(Cds_S) + radical(C=CJO) + radical(Cds_P)"""),
)

species(
    label = '[CH]=C=[C]OC[C]=C[O](25190)',
    structure = SMILES('[CH]=C=[C]OC[C]=C[O]'),
    E0 = (629.363,'kJ/mol'),
    modes = [
        HarmonicOscillator(frequencies=([3120,650,792.5,1650,2750,2850,1437.5,1250,1305,750,350,540,610,2055,3010,987.5,1337.5,450,1655,1670,1700,300,440,180,180,180,180],'cm^-1')),
        HinderedRotor(inertia=(1.12076,'amu*angstrom^2'), symmetry=1, barrier=(25.7685,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(1.12084,'amu*angstrom^2'), symmetry=1, barrier=(25.7704,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(1.12278,'amu*angstrom^2'), symmetry=1, barrier=(25.8149,'kJ/mol'), semiclassical=False),
    ],
    spinMultiplicity = 5,
    opticalIsomers = 1,
    molecularWeight = (108.095,'amu'),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[-0.282264,0.0812221,-9.52014e-05,5.28023e-08,-1.09247e-11,75860.3,33.5095], Tmin=(100,'K'), Tmax=(1334.35,'K')), NASAPolynomial(coeffs=[21.6125,0.0059898,1.60123e-07,-2.32911e-10,2.18005e-14,70871.7,-75.2334], Tmin=(1334.35,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(629.363,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(270.22,'J/(mol*K)'), comment="""Thermo group additivity estimation: group(O2s-Cs(Cds-Cd)) + group(O2s-(Cds-Cd)H) + group(Cs-(Cds-Cds)OsHH) + group(Cds-CdsCsH) + group(Cds-CdsOsH) + group(Cds-CdsOsH) + group(Cds-CdsHH) + group(Cdd-CdsCds) + radical(C=COJ) + radical(C=C=CJ) + radical(Cds_S) + radical(C=CJO)"""),
)

species(
    label = '[CH]=C[CH]O[C]=C=C[O](25191)',
    structure = SMILES('[CH]=C[CH]O[C]=C=C[O]'),
    E0 = (595.08,'kJ/mol'),
    modes = [
        HarmonicOscillator(frequencies=([1685,370,3025,407.5,1350,352.5,2995,3025,975,1000,1300,1375,400,500,1630,1680,3120,650,792.5,1650,540,610,2055,180,180,180,180],'cm^-1')),
        HinderedRotor(inertia=(1.33366,'amu*angstrom^2'), symmetry=1, barrier=(30.6635,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(1.33263,'amu*angstrom^2'), symmetry=1, barrier=(30.6398,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(1.333,'amu*angstrom^2'), symmetry=1, barrier=(30.6483,'kJ/mol'), semiclassical=False),
    ],
    spinMultiplicity = 5,
    opticalIsomers = 1,
    molecularWeight = (108.095,'amu'),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[-0.100944,0.0763449,-8.17859e-05,4.16239e-08,-8.0121e-12,71731.1,33.1078], Tmin=(100,'K'), Tmax=(1383.5,'K')), NASAPolynomial(coeffs=[21.5705,0.00817391,-1.8959e-06,2.46459e-10,-1.45784e-14,66262.3,-76.5719], Tmin=(1383.5,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(595.08,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(270.22,'J/(mol*K)'), comment="""Thermo group additivity estimation: group(O2s-Cs(Cds-Cd)) + group(O2s-(Cds-Cd)H) + group(Cs-(Cds-Cds)OsHH) + group(Cds-CdsCsH) + group(Cds-CdsOsH) + group(Cds-CdsHH) + group(Cds-CdsOsH) + group(Cdd-CdsCds) + radical(Cds_P) + radical(C=CJO) + radical(C=COJ) + radical(C=CCJ(O)C)"""),
)

species(
    label = '[CH]=C=[C]O[CH]C=C[O](25192)',
    structure = SMILES('[CH]=C=[C]O[CH]C=C[O]'),
    E0 = (502.461,'kJ/mol'),
    modes = [
        HarmonicOscillator(frequencies=([3120,650,792.5,1650,3025,407.5,1350,352.5,540,610,2055,2995,3025,975,1000,1300,1375,400,500,1630,1680,1685,370,180,180,180,180],'cm^-1')),
        HinderedRotor(inertia=(1.51234,'amu*angstrom^2'), symmetry=1, barrier=(34.7717,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(1.51252,'amu*angstrom^2'), symmetry=1, barrier=(34.7759,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(1.51586,'amu*angstrom^2'), symmetry=1, barrier=(34.8527,'kJ/mol'), semiclassical=False),
    ],
    spinMultiplicity = 5,
    opticalIsomers = 1,
    molecularWeight = (108.095,'amu'),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[-0.222423,0.0763096,-8.1955e-05,4.19313e-08,-8.01413e-12,60598.3,33.7649], Tmin=(100,'K'), Tmax=(1458.21,'K')), NASAPolynomial(coeffs=[21.2079,0.0073384,-5.29446e-07,-8.51993e-11,1.05479e-14,55431.2,-73.9944], Tmin=(1458.21,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(502.461,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(270.22,'J/(mol*K)'), comment="""Thermo group additivity estimation: group(O2s-Cs(Cds-Cd)) + group(O2s-(Cds-Cd)H) + group(Cs-(Cds-Cds)OsHH) + group(Cds-CdsCsH) + group(Cds-CdsOsH) + group(Cds-CdsOsH) + group(Cds-CdsHH) + group(Cdd-CdsCds) + radical(C=C=CJ) + radical(C=CCJ(O)C) + radical(C=CJO) + radical(C=COJ)"""),
)

species(
    label = '[CH]=C=[C]OC=[C]C[O](25193)',
    structure = SMILES('[CH]=C=[C]OC=[C]C[O]'),
    E0 = (778.999,'kJ/mol'),
    modes = [
        HarmonicOscillator(frequencies=([3120,650,792.5,1650,2750,2850,1437.5,1250,1305,750,350,540,610,2055,3010,987.5,1337.5,450,1655,1670,1700,300,440,180,180,180,1924.69],'cm^-1')),
        HinderedRotor(inertia=(0.401546,'amu*angstrom^2'), symmetry=1, barrier=(9.23233,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.401381,'amu*angstrom^2'), symmetry=1, barrier=(9.22855,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(2.78378,'amu*angstrom^2'), symmetry=1, barrier=(64.0045,'kJ/mol'), semiclassical=False),
    ],
    spinMultiplicity = 5,
    opticalIsomers = 1,
    molecularWeight = (108.095,'amu'),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[0.752292,0.0797457,-0.000129684,1.1399e-07,-3.85139e-11,93801,34.2069], Tmin=(100,'K'), Tmax=(880.552,'K')), NASAPolynomial(coeffs=[7.73549,0.0297697,-1.34563e-05,2.45093e-09,-1.62643e-13,93278.8,5.42366], Tmin=(880.552,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(778.999,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(270.22,'J/(mol*K)'), comment="""Thermo group additivity estimation: group(O2s-(Cds-Cd)(Cds-Cd)) + group(O2s-CsH) + group(Cs-(Cds-Cds)OsHH) + group(Cds-CdsCsH) + group(Cds-CdsOsH) + group(Cds-CdsOsH) + group(Cds-CdsHH) + group(Cdd-CdsCds) + radical(C=C=CJ) + radical(CCOJ) + radical(Cds_S) + radical(C=CJO)"""),
)

species(
    label = '[CH]=[C]CO[CH][C]=C=O(25194)',
    structure = SMILES('[CH]=[C]CO[CH][C]=C=O'),
    E0 = (688.95,'kJ/mol'),
    modes = [
        HarmonicOscillator(frequencies=([1670,1700,300,440,2750,2850,1437.5,1250,1305,750,350,2120,512.5,787.5,3120,650,792.5,1650,3025,407.5,1350,352.5,180,180,180,1355.82],'cm^-1')),
        HinderedRotor(inertia=(0.218733,'amu*angstrom^2'), symmetry=1, barrier=(5.02911,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.218471,'amu*angstrom^2'), symmetry=1, barrier=(5.02307,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.218404,'amu*angstrom^2'), symmetry=1, barrier=(5.02154,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(2.44335,'amu*angstrom^2'), symmetry=1, barrier=(56.1775,'kJ/mol'), semiclassical=False),
    ],
    spinMultiplicity = 5,
    opticalIsomers = 1,
    molecularWeight = (108.095,'amu'),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[0.476462,0.0930627,-0.000177103,1.72176e-07,-6.21783e-11,82973.5,31.1487], Tmin=(100,'K'), Tmax=(874.052,'K')), NASAPolynomial(coeffs=[4.60827,0.037523,-1.89248e-05,3.58033e-09,-2.41382e-13,83650.5,19.7757], Tmin=(874.052,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(688.95,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(266.063,'J/(mol*K)'), comment="""Thermo group additivity estimation: group(O2s-CsCs) + group(Cs-(Cds-Cds)OsHH) + group(Cs-(Cds-Cdd-O2d)OsHH) + group(Cds-CdsCsH) + group(Cds-(Cdd-O2d)CsH) + group(Cds-CdsHH) + radical(C=CCJ(O)C) + radical(Cds_S) + radical(Cds_S) + radical(Cds_P)"""),
)

species(
    label = '[CH]C=CO[CH][C]=C=O(25195)',
    structure = SMILES('[CH]C=CO[CH][C]=C=O'),
    E0 = (516.978,'kJ/mol'),
    modes = [
        HarmonicOscillator(frequencies=([3025,407.5,1350,352.5,2120,512.5,787.5,2995,3025,975,1000,1300,1375,400,500,1630,1680,1685,370,200,800,960,1120,1280,1440,1600],'cm^-1')),
        HinderedRotor(inertia=(0.156089,'amu*angstrom^2'), symmetry=1, barrier=(3.5888,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.156089,'amu*angstrom^2'), symmetry=1, barrier=(3.5888,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.156089,'amu*angstrom^2'), symmetry=1, barrier=(3.5888,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.156089,'amu*angstrom^2'), symmetry=1, barrier=(3.5888,'kJ/mol'), semiclassical=False),
    ],
    spinMultiplicity = 5,
    opticalIsomers = 1,
    molecularWeight = (108.095,'amu'),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[0.202671,0.0823703,-9.44759e-05,5.6224e-08,-1.3251e-11,62315.8,28.7702], Tmin=(100,'K'), Tmax=(1036.92,'K')), NASAPolynomial(coeffs=[15.3837,0.0238084,-9.76067e-06,1.75794e-09,-1.19297e-13,59167.6,-45.0197], Tmin=(1036.92,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(516.978,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(266.063,'J/(mol*K)'), comment="""Thermo group additivity estimation: group(O2s-Cs(Cds-Cd)) + group(Cs-(Cds-Cdd-O2d)OsHH) + group(Cs-(Cds-Cds)HHH) + group(Cds-CdsCsH) + group(Cds-(Cdd-O2d)CsH) + group(Cds-CdsOsH) + radical(Cds_S) + radical(AllylJ2_triplet) + radical(C=CCJ(O)C)"""),
)

species(
    label = '[CH]=C=C[C](C=O)C=O(25196)',
    structure = SMILES('[CH]=C=C[C](C=O)C=O'),
    E0 = (202.26,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
    molecularWeight = (108.095,'amu'),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[1.54055,0.0448725,-1.3526e-05,-1.42212e-08,8.4592e-12,24422.8,26.9687], Tmin=(100,'K'), Tmax=(1017.36,'K')), NASAPolynomial(coeffs=[12.4462,0.0217391,-8.52982e-06,1.58157e-09,-1.11907e-13,21182,-30.8544], Tmin=(1017.36,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(202.26,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(270.22,'J/(mol*K)'), comment="""Thermo group additivity estimation: group(Cs-CsCsCsH) + group(Cds-CdsCsH) + group(Cds-OdCsH) + group(Cds-OdCsH) + group(Cds-CdsHH) + group(Cdd-CdsCds) + radical(C=C=CJ) + radical(C=CCJ(C)C=O)"""),
)

species(
    label = '[CH]=C(C=O)C=C=C[O](22652)',
    structure = SMILES('[CH]=C(C=O)C=C=C[O]'),
    E0 = (300.313,'kJ/mol'),
    modes = [
        HarmonicOscillator(frequencies=([3120,650,792.5,1650,2782.5,750,1395,475,1775,1000,540,610,2055,350,440,435,1725,2995,3025,975,1000,1300,1375,400,500,1630,1680,180],'cm^-1')),
        HinderedRotor(inertia=(1.26441,'amu*angstrom^2'), symmetry=1, barrier=(29.0713,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(1.26476,'amu*angstrom^2'), symmetry=1, barrier=(29.0794,'kJ/mol'), semiclassical=False),
    ],
    spinMultiplicity = 3,
    opticalIsomers = 1,
    molecularWeight = (108.095,'amu'),
    collisionModel = TransportData(shapeIndex=2, epsilon=(4313.19,'J/mol'), sigma=(6.61299,'angstroms'), dipoleMoment=(0,'C*m'), polarizability=(0,'angstroms^3'), rotrelaxcollnum=0, comment="""Epsilon & sigma estimated with Tc=673.71 K, Pc=33.84 bar (from Joback method)"""),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[0.0723359,0.0739509,-6.32921e-05,1.35797e-08,3.91554e-12,36271.8,25.8627], Tmin=(100,'K'), Tmax=(977.162,'K')), NASAPolynomial(coeffs=[23.1468,0.00679682,-2.11532e-06,4.34057e-10,-3.63074e-14,30458.9,-91.5945], Tmin=(977.162,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(300.313,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(274.378,'J/(mol*K)'), comment="""Thermo group additivity estimation: group(O2s-(Cds-Cd)H) + group(Cds-Cds(Cds-O2d)(Cds-Cds)) + group(Cds-Cds(Cds-Cds)H) + group(Cds-O2d(Cds-Cds)H) + group(Cds-CdsHH) + group(Cds-CdsOsH) + group(Cdd-CdsCds) + radical(Cds_P) + radical(C=COJ)"""),
)

species(
    label = '[CH]C1=CO[CH]C1=C[O](25197)',
    structure = SMILES('[CH]C1=CO[CH]C1=C[O]'),
    E0 = (336.041,'kJ/mol'),
    modes = [
        HarmonicOscillator(frequencies=([2750,3150,900,1100,3010,987.5,1337.5,450,1655,300,800,800,800,800,800,800,800,800,800,1600,1600,1600,1600,1600,1600,1600,1600,1600,1600],'cm^-1')),
        HinderedRotor(inertia=(0.156089,'amu*angstrom^2'), symmetry=1, barrier=(3.5888,'kJ/mol'), semiclassical=False),
    ],
    spinMultiplicity = 5,
    opticalIsomers = 1,
    molecularWeight = (108.095,'amu'),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[1.04554,0.0302925,0.000101154,-1.78015e-07,7.7807e-11,40555.1,23.4574], Tmin=(100,'K'), Tmax=(914.954,'K')), NASAPolynomial(coeffs=[29.1851,-0.00244477,6.81125e-06,-1.42634e-09,8.83808e-14,31626.8,-130.45], Tmin=(914.954,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(336.041,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(278.535,'J/(mol*K)'), comment="""Thermo group additivity estimation: group(O2s-Cs(Cds-Cd)) + group(O2s-(Cds-Cd)H) + group(Cs-(Cds-Cds)OsHH) + group(Cs-(Cds-Cds)HHH) + group(Cds-Cds(Cds-Cds)Cs) + group(Cds-Cds(Cds-Cds)Cs) + group(Cds-CdsOsH) + group(Cds-CdsOsH) + ring(Cyclopentane) + radical(C=CCJ(O)C) + radical(AllylJ2_triplet) + radical(C=COJ)"""),
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
    label = '[CH]=C=COC=C=[CH](22632)',
    structure = SMILES('[CH]=C=COC=C=[CH]'),
    E0 = (559.925,'kJ/mol'),
    modes = [
        HarmonicOscillator(frequencies=([3115,3125,620,680,785,800,1600,1700,2995,3025,975,1000,1300,1375,400,500,1630,1680,540,563.333,586.667,610,1970,2140,180],'cm^-1')),
        HinderedRotor(inertia=(1.58688,'amu*angstrom^2'), symmetry=1, barrier=(36.4854,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(1.58588,'amu*angstrom^2'), symmetry=1, barrier=(36.4626,'kJ/mol'), semiclassical=False),
    ],
    spinMultiplicity = 3,
    opticalIsomers = 1,
    molecularWeight = (92.0954,'amu'),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[0.16531,0.069971,-7.94956e-05,4.33001e-08,-8.70502e-12,67493.8,27.0755], Tmin=(100,'K'), Tmax=(1426.9,'K')), NASAPolynomial(coeffs=[18.4934,0.00603881,9.08569e-07,-4.31438e-10,3.72195e-14,63541.3,-63.3845], Tmin=(1426.9,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(559.925,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(249.434,'J/(mol*K)'), comment="""Thermo group additivity estimation: group(O2s-(Cds-Cd)(Cds-Cd)) + group(Cds-CdsOsH) + group(Cds-CdsOsH) + group(Cds-CdsHH) + group(Cds-CdsHH) + group(Cdd-CdsCds) + group(Cdd-CdsCds) + radical(C=C=CJ) + radical(C=C=CJ)"""),
)

species(
    label = '[C]#C[CH]OC=C=C[O](25198)',
    structure = SMILES('[C]#C[CH]OC=C=C[O]'),
    E0 = (694.322,'kJ/mol'),
    modes = [
        HarmonicOscillator(frequencies=([3025,407.5,1350,352.5,2995,3025,975,1000,1300,1375,400,500,1630,1680,2175,525,540,610,2055,180,180,180,180,180],'cm^-1')),
        HinderedRotor(inertia=(1.42642,'amu*angstrom^2'), symmetry=1, barrier=(32.7962,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(1.42767,'amu*angstrom^2'), symmetry=1, barrier=(32.825,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(1.42802,'amu*angstrom^2'), symmetry=1, barrier=(32.833,'kJ/mol'), semiclassical=False),
    ],
    spinMultiplicity = 4,
    opticalIsomers = 1,
    molecularWeight = (107.087,'amu'),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[-1.23665,0.0944963,-0.000121629,6.93727e-08,-1.42459e-11,83714.1,30.5289], Tmin=(100,'K'), Tmax=(1421.37,'K')), NASAPolynomial(coeffs=[26.5869,-0.00546332,6.71781e-06,-1.54658e-09,1.13446e-13,77992.4,-105.791], Tmin=(1421.37,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(694.322,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(245.277,'J/(mol*K)'), comment="""Thermo group additivity estimation: group(O2s-Cs(Cds-Cd)) + group(O2s-(Cds-Cd)H) + group(Cs-CtOsHH) + group(Cds-CdsOsH) + group(Cds-CdsOsH) + group(Cdd-CdsCds) + group(Ct-CtCs) + group(Ct-CtH) + radical(C=COJ) + radical(Acetyl) + radical(CCsJOC(O))"""),
)

species(
    label = 'C#CC1OC1[C]=C[O](25199)',
    structure = SMILES('C#CC1OC1[C]=C[O]'),
    E0 = (403.94,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
    molecularWeight = (108.095,'amu'),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[-1.36598,0.0778268,-8.14651e-05,4.03513e-08,-7.09871e-12,48810.9,33.8742], Tmin=(100,'K'), Tmax=(1764.1,'K')), NASAPolynomial(coeffs=[17.667,0.00601892,3.95478e-06,-1.13625e-09,8.54507e-14,46553.9,-56.1169], Tmin=(1764.1,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(403.94,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(274.378,'J/(mol*K)'), comment="""Thermo group additivity estimation: group(O2s-CsCs) + group(O2s-(Cds-Cd)H) + group(Cs-(Cds-Cds)CsOsH) + group(Cs-(Cds-Cds)CsOsH) + group(Cds-CdsCsH) + group(Cds-CdsOsH) + group(Ct-CtCs) + group(Ct-CtH) + ring(Ethylene_oxide) + radical(Cds_S) + radical(C=COJ)"""),
)

species(
    label = 'C#CCO[C]=C=C[O](25200)',
    structure = SMILES('C#CCO[C]=C=C[O]'),
    E0 = (402.995,'kJ/mol'),
    modes = [
        HarmonicOscillator(frequencies=([750,770,3400,2100,540,610,2055,2175,525,2750,2850,1437.5,1250,1305,750,350,1685,370,3010,987.5,1337.5,450,1655,180,180,180,180],'cm^-1')),
        HinderedRotor(inertia=(1.42489,'amu*angstrom^2'), symmetry=1, barrier=(32.7611,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(1.42624,'amu*angstrom^2'), symmetry=1, barrier=(32.792,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(1.42493,'amu*angstrom^2'), symmetry=1, barrier=(32.7619,'kJ/mol'), semiclassical=False),
    ],
    spinMultiplicity = 3,
    opticalIsomers = 1,
    molecularWeight = (108.095,'amu'),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[-0.571518,0.0818063,-9.28244e-05,4.91395e-08,-9.57677e-12,48649.8,31.0679], Tmin=(100,'K'), Tmax=(1463.22,'K')), NASAPolynomial(coeffs=[22.7658,0.00382449,1.65894e-06,-5.33877e-10,4.21921e-14,43338.8,-85.2153], Tmin=(1463.22,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(402.995,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(270.22,'J/(mol*K)'), comment="""Thermo group additivity estimation: group(O2s-Cs(Cds-Cd)) + group(O2s-(Cds-Cd)H) + group(Cs-CtOsHH) + group(Cds-CdsOsH) + group(Cds-CdsOsH) + group(Cdd-CdsCds) + group(Ct-CtCs) + group(Ct-CtH) + radical(C=COJ) + radical(C=CJO)"""),
)

species(
    label = '[C]#CCOC=C=C[O](25201)',
    structure = SMILES('[C]#CCOC=C=C[O]'),
    E0 = (500.394,'kJ/mol'),
    modes = [
        HarmonicOscillator(frequencies=([2750,2850,1437.5,1250,1305,750,350,2995,3025,975,1000,1300,1375,400,500,1630,1680,2175,525,540,610,2055,180,180,180,180,180],'cm^-1')),
        HinderedRotor(inertia=(1.49228,'amu*angstrom^2'), symmetry=1, barrier=(34.3105,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(1.49982,'amu*angstrom^2'), symmetry=1, barrier=(34.4839,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(1.49893,'amu*angstrom^2'), symmetry=1, barrier=(34.4632,'kJ/mol'), semiclassical=False),
    ],
    spinMultiplicity = 3,
    opticalIsomers = 1,
    molecularWeight = (108.095,'amu'),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[-1.05278,0.087089,-0.00010093,5.36982e-08,-1.03579e-11,60386.2,30.4072], Tmin=(100,'K'), Tmax=(1519.6,'K')), NASAPolynomial(coeffs=[24.0539,0.000906024,3.97738e-06,-1.02816e-09,7.71654e-14,55075.9,-93.5899], Tmin=(1519.6,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(500.394,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(270.22,'J/(mol*K)'), comment="""Thermo group additivity estimation: group(O2s-Cs(Cds-Cd)) + group(O2s-(Cds-Cd)H) + group(Cs-CtOsHH) + group(Cds-CdsOsH) + group(Cds-CdsOsH) + group(Cdd-CdsCds) + group(Ct-CtCs) + group(Ct-CtH) + radical(Acetyl) + radical(C=COJ)"""),
)

species(
    label = 'C#CCO[CH][C]=C=O(25202)',
    structure = SMILES('C#CCO[CH][C]=C=O'),
    E0 = (369.962,'kJ/mol'),
    modes = [
        HarmonicOscillator(frequencies=([750,770,3400,2100,3025,407.5,1350,352.5,2750,2850,1437.5,1250,1305,750,350,2120,512.5,787.5,1685,370,2175,525,180,180,180,800.392],'cm^-1')),
        HinderedRotor(inertia=(0.138462,'amu*angstrom^2'), symmetry=1, barrier=(3.18352,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.145949,'amu*angstrom^2'), symmetry=1, barrier=(3.35565,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(2.52099,'amu*angstrom^2'), symmetry=1, barrier=(57.9625,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(2.51881,'amu*angstrom^2'), symmetry=1, barrier=(57.9124,'kJ/mol'), semiclassical=False),
    ],
    spinMultiplicity = 3,
    opticalIsomers = 1,
    molecularWeight = (108.095,'amu'),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[0.534209,0.0880735,-0.000155326,1.42879e-07,-4.94027e-11,44609.7,27.6829], Tmin=(100,'K'), Tmax=(890.567,'K')), NASAPolynomial(coeffs=[6.61208,0.0328344,-1.52256e-05,2.773e-09,-1.82685e-13,44635.1,5.28577], Tmin=(890.567,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(369.962,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(266.063,'J/(mol*K)'), comment="""Thermo group additivity estimation: group(O2s-CsCs) + group(Cs-(Cds-Cdd-O2d)OsHH) + group(Cs-CtOsHH) + group(Cds-(Cdd-O2d)CsH) + group(Ct-CtCs) + group(Ct-CtH) + radical(Cds_S) + radical(C=CCJ(O)C)"""),
)

species(
    label = '[C]#C[CH]OC=C=CO(25203)',
    structure = SMILES('[C]#C[CH]OC=C=CO'),
    E0 = (552.859,'kJ/mol'),
    modes = [
        HarmonicOscillator(frequencies=([3025,407.5,1350,352.5,540,610,2055,3615,1277.5,1000,2175,525,2995,3025,975,1000,1300,1375,400,500,1630,1680,180,180,180,180],'cm^-1')),
        HinderedRotor(inertia=(1.28083,'amu*angstrom^2'), symmetry=1, barrier=(29.4489,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(1.28127,'amu*angstrom^2'), symmetry=1, barrier=(29.4589,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(1.28346,'amu*angstrom^2'), symmetry=1, barrier=(29.5092,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(1.27817,'amu*angstrom^2'), symmetry=1, barrier=(29.3876,'kJ/mol'), semiclassical=False),
    ],
    spinMultiplicity = 3,
    opticalIsomers = 1,
    molecularWeight = (108.095,'amu'),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[-1.8928,0.104407,-0.000135603,7.73243e-08,-1.57613e-11,66727.6,31.4701], Tmin=(100,'K'), Tmax=(1448.41,'K')), NASAPolynomial(coeffs=[29.1295,-0.00818348,8.87474e-06,-2.00587e-09,1.45998e-13,60564.5,-119.94], Tmin=(1448.41,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(552.859,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(266.063,'J/(mol*K)'), comment="""Thermo group additivity estimation: group(O2s-Cs(Cds-Cd)) + group(O2s-(Cds-Cd)H) + group(Cs-CtOsHH) + group(Cds-CdsOsH) + group(Cds-CdsOsH) + group(Cdd-CdsCds) + group(Ct-CtCs) + group(Ct-CtH) + radical(CCsJOC(O)) + radical(Acetyl)"""),
)

species(
    label = 'C#CC1O[CH]C1=C[O](25204)',
    structure = SMILES('C#CC1O[CH]C1=C[O]'),
    E0 = (277.048,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
    molecularWeight = (108.095,'amu'),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[1.36797,0.0397035,2.49713e-05,-6.9879e-08,3.19577e-11,33432.6,25.2243], Tmin=(100,'K'), Tmax=(946.772,'K')), NASAPolynomial(coeffs=[18.7539,0.0117725,-2.89929e-06,5.30785e-10,-4.43429e-14,28100.2,-68.4765], Tmin=(946.772,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(277.048,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(278.535,'J/(mol*K)'), comment="""Thermo group additivity estimation: group(O2s-CsCs) + group(O2s-(Cds-Cd)H) + group(Cs-(Cds-Cds)(Cds-Cds)OsH) + group(Cs-(Cds-Cds)OsHH) + group(Cds-CdsCsCs) + group(Cds-CdsOsH) + group(Ct-CtCs) + group(Ct-CtH) + ring(Cyclobutane) + radical(C=COJ) + radical(C=CCJ(O)C)"""),
)

species(
    label = 'C#CC1OC=[C]C1[O](25205)',
    structure = SMILES('C#CC1OC=[C]C1[O]'),
    E0 = (400.065,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
    molecularWeight = (108.095,'amu'),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[1.17806,0.0381003,4.63305e-05,-1.12427e-07,5.41093e-11,48241,21.7781], Tmin=(100,'K'), Tmax=(891.765,'K')), NASAPolynomial(coeffs=[25.1196,-0.00415455,7.84559e-06,-1.75103e-09,1.20688e-13,41381,-105.504], Tmin=(891.765,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(400.065,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(278.535,'J/(mol*K)'), comment="""Thermo group additivity estimation: group(O2s-Cs(Cds-Cd)) + group(O2s-CsH) + group(Cs-(Cds-Cds)CsOsH) + group(Cs-(Cds-Cds)CsOsH) + group(Cds-CdsCsH) + group(Cds-CdsOsH) + group(Ct-CtCs) + group(Ct-CtH) + ring(2,3-Dihydrofuran) + radical(CC(C)OJ) + radical(Cds_S)"""),
)

species(
    label = '[C]1[CH]OC=C=COC=1(25206)',
    structure = SMILES('[C]1[CH]OC=C=COC=1'),
    E0 = (345.541,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
    molecularWeight = (108.095,'amu'),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[1.85744,0.0164796,0.000108095,-1.66475e-07,6.87744e-11,41664.4,18.1026], Tmin=(100,'K'), Tmax=(932.13,'K')), NASAPolynomial(coeffs=[23.0953,0.00381511,2.19613e-06,-4.19462e-10,1.49945e-14,34296.1,-101.152], Tmin=(932.13,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(345.541,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(282.692,'J/(mol*K)'), comment="""Thermo group additivity estimation: group(O2s-Cs(Cds-Cd)) + group(O2s-(Cds-Cd)(Cds-Cd)) + group(Cs-(Cds-Cds)OsHH) + group(Cds-CdsCsH) + group(Cds-CdsOsH) + group(Cds-CdsOsH) + group(Cds-CdsOsH) + group(Cdd-CdsCds) + ring(Cyclooctane) + radical(Cds_S) + radical(C=CCJ(O)C)"""),
)

species(
    label = 'C#CCOC=C=C=O(25207)',
    structure = SMILES('C#CCOC=C=C=O'),
    E0 = (207.347,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
    molecularWeight = (108.095,'amu'),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[1.27083,0.0473408,-1.54107e-05,-2.98768e-08,1.98487e-11,25048.3,8.53122], Tmin=(100,'K'), Tmax=(904.312,'K')), NASAPolynomial(coeffs=[19.335,0.00206475,2.25426e-06,-5.57761e-10,3.77864e-14,20365.3,-84.6292], Tmin=(904.312,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(207.347,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(270.22,'J/(mol*K)'), comment="""Thermo group additivity estimation: group(O2s-Cs(Cds-Cd)) + group(Cs-CtOsHH) + group(Cds-CdsOsH) + group(Ct-CtCs) + group(Cdd-CdsCds) + group(Ct-CtH)"""),
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
    label = '[CH]OC=C=C[O](23145)',
    structure = SMILES('[CH]OC=C=C[O]'),
    E0 = (386.278,'kJ/mol'),
    modes = [
        HarmonicOscillator(frequencies=([2995,3025,975,1000,1300,1375,400,500,1630,1680,540,610,2055,180,180,180,180,180,180],'cm^-1')),
        HinderedRotor(inertia=(1.0316,'amu*angstrom^2'), symmetry=1, barrier=(23.7185,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(1.03202,'amu*angstrom^2'), symmetry=1, barrier=(23.7281,'kJ/mol'), semiclassical=False),
    ],
    spinMultiplicity = 4,
    opticalIsomers = 1,
    molecularWeight = (83.0654,'amu'),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[0.947942,0.0502343,-1.58912e-05,-3.98424e-08,2.5868e-11,46584.2,19.7696], Tmin=(100,'K'), Tmax=(908.33,'K')), NASAPolynomial(coeffs=[24.7588,-0.00908885,6.88222e-06,-1.37011e-09,9.0156e-14,40380.2,-103.154], Tmin=(908.33,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(386.278,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(199.547,'J/(mol*K)'), comment="""Thermo group additivity estimation: group(O2s-Cs(Cds-Cd)) + group(O2s-(Cds-Cd)H) + group(Cs-OsHHH) + group(Cds-CdsOsH) + group(Cds-CdsOsH) + group(Cdd-CdsCds) + radical(C=COJ) + radical(CH2_triplet)"""),
)

species(
    label = '[CH]=C=COC#CC=O(25208)',
    structure = SMILES('[CH]=C=COC#CC=O'),
    E0 = (344.235,'kJ/mol'),
    modes = [
        HarmonicOscillator(frequencies=([3120,650,792.5,1650,2782.5,750,1395,475,1775,1000,3010,987.5,1337.5,450,1655,2100,2250,500,550,540,610,2055,180,180],'cm^-1')),
        HinderedRotor(inertia=(1.45123,'amu*angstrom^2'), symmetry=1, barrier=(33.3665,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(1.45073,'amu*angstrom^2'), symmetry=1, barrier=(33.3551,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(1.45319,'amu*angstrom^2'), symmetry=1, barrier=(33.4116,'kJ/mol'), semiclassical=False),
    ],
    spinMultiplicity = 2,
    opticalIsomers = 1,
    molecularWeight = (107.087,'amu'),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[0.712453,0.0694977,-8.06658e-05,4.65682e-08,-1.04659e-11,41522.7,26.6106], Tmin=(100,'K'), Tmax=(1094.1,'K')), NASAPolynomial(coeffs=[15.6007,0.0150657,-6.0389e-06,1.09518e-09,-7.52188e-14,38264.9,-46.5553], Tmin=(1094.1,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(344.235,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(245.277,'J/(mol*K)'), comment="""Thermo group additivity estimation: group(O2s-(Cds-Cd)(Cds-Cd)) + group(Cds-CdsOsH) + group(Cds-CdsHH) + group(Cds-O2d(Cds-Cds)H) + group(Cdd-CdsCds) + group(Ct-CtCs) + group(Ct-CtOs) + radical(C=C=CJ)"""),
)

species(
    label = '[CH]=C=CO[C]=CC=O(25209)',
    structure = SMILES('[CH]=C=CO[C]=CC=O'),
    E0 = (381.125,'kJ/mol'),
    modes = [
        HarmonicOscillator(frequencies=([3120,650,792.5,1650,2782.5,750,1395,475,1775,1000,540,610,2055,2995,3025,975,1000,1300,1375,400,500,1630,1680,1685,370,180,180],'cm^-1')),
        HinderedRotor(inertia=(0.951645,'amu*angstrom^2'), symmetry=1, barrier=(21.8802,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.950827,'amu*angstrom^2'), symmetry=1, barrier=(21.8614,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.951026,'amu*angstrom^2'), symmetry=1, barrier=(21.8659,'kJ/mol'), semiclassical=False),
    ],
    spinMultiplicity = 3,
    opticalIsomers = 1,
    molecularWeight = (108.095,'amu'),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[0.909155,0.0712494,-8.61062e-05,5.5043e-08,-1.41775e-11,45947.4,29.8374], Tmin=(100,'K'), Tmax=(941.726,'K')), NASAPolynomial(coeffs=[11.9343,0.0244189,-1.15118e-05,2.23482e-09,-1.58139e-13,43870.9,-22.6901], Tmin=(941.726,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(381.125,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(270.22,'J/(mol*K)'), comment="""Thermo group additivity estimation: group(O2s-(Cds-Cd)(Cds-Cd)) + group(Cd-Cd(CO)H) + group(Cds-CdsOsH) + group(Cds-CdsOsH) + group(Cds-O2d(Cds-Cds)H) + group(Cds-CdsHH) + group(Cdd-CdsCds) + radical(C=C=CJ) + radical(C=CJO)"""),
)

species(
    label = '[CH]=C=CO[CH]C=C=O(25210)',
    structure = SMILES('[CH]=C=CO[CH]C=C=O'),
    E0 = (239.532,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
    molecularWeight = (108.095,'amu'),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[0.159912,0.0800176,-9.24762e-05,4.8897e-08,-8.61409e-12,28951.5,27.0196], Tmin=(100,'K'), Tmax=(901.269,'K')), NASAPolynomial(coeffs=[18.6015,0.0123482,-3.44889e-06,4.97776e-10,-3.03345e-14,25051.6,-63.2278], Tmin=(901.269,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(239.532,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(270.22,'J/(mol*K)'), comment="""Thermo group additivity estimation: group(O2s-Cs(Cds-Cd)) + group(Cs-(Cds-Cdd-O2d)OsHH) + group(Cds-(Cdd-O2d)CsH) + group(Cds-CdsOsH) + group(Cds-CdsHH) + group(Cdd-CdsCds) + radical(C=C=CJ) + radical(C=CCJ(O)C)"""),
)

species(
    label = '[CH]=C=[C]OC=CC=O(25211)',
    structure = SMILES('[CH]=C=[C]OC=CC=O'),
    E0 = (381.125,'kJ/mol'),
    modes = [
        HarmonicOscillator(frequencies=([3120,650,792.5,1650,2782.5,750,1395,475,1775,1000,540,610,2055,2995,3025,975,1000,1300,1375,400,500,1630,1680,1685,370,180,180],'cm^-1')),
        HinderedRotor(inertia=(0.951645,'amu*angstrom^2'), symmetry=1, barrier=(21.8802,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.950827,'amu*angstrom^2'), symmetry=1, barrier=(21.8614,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.951026,'amu*angstrom^2'), symmetry=1, barrier=(21.8659,'kJ/mol'), semiclassical=False),
    ],
    spinMultiplicity = 3,
    opticalIsomers = 1,
    molecularWeight = (108.095,'amu'),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[0.909155,0.0712494,-8.61062e-05,5.5043e-08,-1.41775e-11,45947.4,29.8374], Tmin=(100,'K'), Tmax=(941.726,'K')), NASAPolynomial(coeffs=[11.9343,0.0244189,-1.15118e-05,2.23482e-09,-1.58139e-13,43870.9,-22.6901], Tmin=(941.726,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(381.125,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(270.22,'J/(mol*K)'), comment="""Thermo group additivity estimation: group(O2s-(Cds-Cd)(Cds-Cd)) + group(Cd-Cd(CO)H) + group(Cds-CdsOsH) + group(Cds-CdsOsH) + group(Cds-O2d(Cds-Cds)H) + group(Cds-CdsHH) + group(Cdd-CdsCds) + radical(C=C=CJ) + radical(C=CJO)"""),
)

species(
    label = '[CH]C1=COC=C1C=O(25071)',
    structure = SMILES('[CH]C1=COC=C1C=O'),
    E0 = (171.942,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
    molecularWeight = (108.095,'amu'),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[1.13394,0.0471244,6.97249e-06,-4.45222e-08,2.00883e-11,20797,23.2298], Tmin=(100,'K'), Tmax=(1014.54,'K')), NASAPolynomial(coeffs=[16.9762,0.0214031,-9.31781e-06,1.87622e-09,-1.40677e-13,15691.7,-62.7468], Tmin=(1014.54,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(171.942,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(274.378,'J/(mol*K)'), comment="""Thermo group additivity estimation: group(O2s-(Cds-Cd)(Cds-Cd)) + group(Cs-(Cds-Cds)HHH) + group(Cds-Cds(Cds-Cds)Cs) + group(Cds-Cds(Cds-O2d)(Cds-Cds)) + group(Cds-CdsOsH) + group(Cds-CdsOsH) + group(Cds-O2d(Cds-Cds)H) + ring(Furan) + radical(AllylJ2_triplet)"""),
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
    label = '[C]=COC=C=[CH](23405)',
    structure = SMILES('[C]=COC=C=[CH]'),
    E0 = (822.971,'kJ/mol'),
    modes = [
        HarmonicOscillator(frequencies=([3120,650,792.5,1650,2995,3025,975,1000,1300,1375,400,500,1630,1680,540,610,2055,180,180],'cm^-1')),
        HinderedRotor(inertia=(1.39111,'amu*angstrom^2'), symmetry=1, barrier=(31.9845,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(1.38597,'amu*angstrom^2'), symmetry=1, barrier=(31.8662,'kJ/mol'), semiclassical=False),
    ],
    spinMultiplicity = 4,
    opticalIsomers = 1,
    molecularWeight = (79.0767,'amu'),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[1.0276,0.0565025,-6.62699e-05,3.70334e-08,-7.71698e-12,99095.2,23.5681], Tmin=(100,'K'), Tmax=(1332.57,'K')), NASAPolynomial(coeffs=[15.9142,0.00468798,7.95695e-08,-1.7505e-10,1.68034e-14,95760.7,-50.1504], Tmin=(1332.57,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(822.971,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(199.547,'J/(mol*K)'), comment="""Thermo group additivity estimation: group(O2s-(Cds-Cd)(Cds-Cd)) + group(Cds-CdsOsH) + group(Cds-CdsOsH) + group(Cds-CdsHH) + group(Cds-CdsHH) + group(Cdd-CdsCds) + radical(C=C=CJ) + radical(CdCdJ2_triplet)"""),
)

species(
    label = '[C]#C[CH]OC=CC=O(25212)',
    structure = SMILES('[C]#C[CH]OC=CC=O'),
    E0 = (497.585,'kJ/mol'),
    modes = [
        HarmonicOscillator(frequencies=([3025,407.5,1350,352.5,2995,3025,975,1000,1300,1375,400,500,1630,1680,2175,525,2782.5,750,1395,475,1775,1000,180,180,180,180],'cm^-1')),
        HinderedRotor(inertia=(1.04104,'amu*angstrom^2'), symmetry=1, barrier=(23.9355,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(1.04086,'amu*angstrom^2'), symmetry=1, barrier=(23.9314,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(1.0401,'amu*angstrom^2'), symmetry=1, barrier=(23.9138,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(1.03959,'amu*angstrom^2'), symmetry=1, barrier=(23.9023,'kJ/mol'), semiclassical=False),
    ],
    spinMultiplicity = 3,
    opticalIsomers = 1,
    molecularWeight = (108.095,'amu'),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[-0.187281,0.0877514,-0.000113805,7.06438e-08,-1.67174e-11,60000.2,27.7236], Tmin=(100,'K'), Tmax=(1050.85,'K')), NASAPolynomial(coeffs=[20.161,0.0102979,-3.24873e-06,5.07021e-10,-3.19288e-14,55723.6,-71.4549], Tmin=(1050.85,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(497.585,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(266.063,'J/(mol*K)'), comment="""Thermo group additivity estimation: group(O2s-Cs(Cds-Cd)) + group(Cs-CtOsHH) + group(Cd-Cd(CO)H) + group(Cds-CdsOsH) + group(Cds-O2d(Cds-Cds)H) + group(Ct-CtCs) + group(Ct-CtH) + radical(Acetyl) + radical(CCsJOC(O))"""),
)

species(
    label = 'C#CC1OC=C1C=O(25213)',
    structure = SMILES('C#CC1OC=C1C=O'),
    E0 = (97.8924,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
    molecularWeight = (108.095,'amu'),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[0.799706,0.0503567,6.64446e-06,-6.01092e-08,3.04223e-11,11907.3,22.061], Tmin=(100,'K'), Tmax=(949.012,'K')), NASAPolynomial(coeffs=[23.4133,0.0048876,-2.72508e-07,9.51149e-11,-1.717e-14,5370.54,-97.6792], Tmin=(949.012,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(97.8924,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(274.378,'J/(mol*K)'), comment="""Thermo group additivity estimation: group(O2s-Cs(Cds-Cd)) + group(Cs-(Cds-Cds)(Cds-Cds)OsH) + group(Cd-CdCs(CO)) + group(Cds-CdsOsH) + group(Cds-O2d(Cds-Cds)H) + group(Ct-CtCs) + group(Ct-CtH) + ring(Cyclobutene)"""),
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
    E0 = (338.118,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS2',
    E0 = (524.68,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS3',
    E0 = (493.202,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS4',
    E0 = (596.941,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS5',
    E0 = (614.085,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS6',
    E0 = (532.735,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS7',
    E0 = (751.561,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS8',
    E0 = (530.818,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS9',
    E0 = (555.941,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS10',
    E0 = (489.639,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS11',
    E0 = (539.749,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS12',
    E0 = (539.749,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS13',
    E0 = (789.667,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS14',
    E0 = (789.667,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS15',
    E0 = (689.179,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS16',
    E0 = (525.317,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS17',
    E0 = (488.116,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS18',
    E0 = (463.776,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS19',
    E0 = (350.511,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS20',
    E0 = (487.062,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS21',
    E0 = (400.288,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS22',
    E0 = (363.091,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS23',
    E0 = (635.067,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS24',
    E0 = (525.322,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS25',
    E0 = (810.951,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS26',
    E0 = (722.457,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS27',
    E0 = (634.305,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS28',
    E0 = (541.686,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS29',
    E0 = (818.224,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS30',
    E0 = (713.923,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS31',
    E0 = (541.951,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS32',
    E0 = (651.918,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS33',
    E0 = (651.918,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS34',
    E0 = (486.353,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS35',
    E0 = (966.806,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS36',
    E0 = (906.126,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS37',
    E0 = (446.274,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS38',
    E0 = (376.762,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS39',
    E0 = (553.409,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS40',
    E0 = (647.671,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS41',
    E0 = (432.439,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS42',
    E0 = (677.295,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS43',
    E0 = (464.757,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS44',
    E0 = (481.086,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS45',
    E0 = (403.538,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS46',
    E0 = (377.343,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS47',
    E0 = (972.53,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS48',
    E0 = (571.149,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS49',
    E0 = (577.585,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS50',
    E0 = (500.897,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS51',
    E0 = (775.43,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS52',
    E0 = (376.506,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS53',
    E0 = (890.648,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS54',
    E0 = (685.939,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS55',
    E0 = (345.649,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

reaction(
    label = 'reaction1',
    reactants = ['[CH]=C=COC=C=C[O](22649)'],
    products = ['C#CC=O(21959)', 'C#CC=O(21959)'],
    transitionState = 'TS1',
    kinetics = Arrhenius(A=(5e+12,'s^-1'), n=0, Ea=(0,'kJ/mol'), T0=(1,'K'), Tmin=(300,'K'), Tmax=(1500,'K'), comment="""Exact match found for rate rule [RJJ]
Euclidian distance = 0
family: 1,4_Linear_birad_scission"""),
)

reaction(
    label = 'reaction2',
    reactants = ['[CH]=C=COC=C=C[O](22649)'],
    products = ['[O]C=[C]C1C=C=CO1(25116)'],
    transitionState = 'TS2',
    kinetics = Arrhenius(A=(7.4226e+09,'s^-1'), n=0.3735, Ea=(186.562,'kJ/mol'), T0=(1,'K'), comment="""Estimated using an average for rate rule [R6;doublebond_intra;radadd_intra_cdsingleH]
Euclidian distance = 0
family: Intra_R_Add_Exocyclic
Ea raised from 183.6 to 186.6 kJ/mol to match endothermicity of reaction."""),
)

reaction(
    label = 'reaction3',
    reactants = ['[CH]=C=COC=C=C[O](22649)'],
    products = ['[CH]=[C]C1OC=C=CO1(25171)'],
    transitionState = 'TS3',
    kinetics = Arrhenius(A=(3.04811e+11,'s^-1'), n=0.222, Ea=(155.084,'kJ/mol'), T0=(1,'K'), comment="""Estimated using average of templates [R7;doublebond_intra;radadd_intra_O] + [R7_SMMS;doublebond_intra;radadd_intra] for rate rule [R7_SMMS;doublebond_intra;radadd_intra_O]
Euclidian distance = 1.0
family: Intra_R_Add_Exocyclic"""),
)

reaction(
    label = 'reaction4',
    reactants = ['H(8)', '[CH]=C=COC=C=C=O(25172)'],
    products = ['[CH]=C=COC=C=C[O](22649)'],
    transitionState = 'TS4',
    kinetics = Arrhenius(A=(7048.5,'m^3/(mol*s)'), n=1.43107, Ea=(2.92155,'kJ/mol'), T0=(1,'K'), Tmin=(303.03,'K'), Tmax=(2000,'K'), comment="""From training reaction 188 used for Ck_O;HJ
Exact match found for rate rule [Ck_O;HJ]
Euclidian distance = 0
family: R_Addition_MultipleBond"""),
)

reaction(
    label = 'reaction5',
    reactants = ['[CH]=C=COC=C=[C]O(25173)'],
    products = ['[CH]=C=COC=C=C[O](22649)'],
    transitionState = 'TS5',
    kinetics = Arrhenius(A=(3.9024e+09,'s^-1'), n=1.09833, Ea=(177.686,'kJ/mol'), T0=(1,'K'), comment="""Estimated using template [R2H_S;Cd_rad_out_double;XH_out] for rate rule [R2H_S;Cd_rad_out_double;O_H_out]
Euclidian distance = 1.0
family: intra_H_migration"""),
)

reaction(
    label = 'reaction6',
    reactants = ['[CH2]C#COC=C=C[O](25174)'],
    products = ['[CH]=C=COC=C=C[O](22649)'],
    transitionState = 'TS6',
    kinetics = Arrhenius(A=(1.65028e+09,'s^-1'), n=1.32317, Ea=(166.12,'kJ/mol'), T0=(1,'K'), comment="""Estimated using average of templates [RnH;Cd_rad_out_singleNd;Cd_H_out_singleH] + [R3H;Cd_rad_out_singleNd;XH_out] for rate rule [R3H;Cd_rad_out_singleNd;Cd_H_out_singleH]
Euclidian distance = 2.0
Multiplied by reaction path degeneracy 2.0
family: intra_H_migration"""),
)

reaction(
    label = 'reaction7',
    reactants = ['[CH]=C=CO[C]=C=CO(25175)'],
    products = ['[CH]=C=COC=C=C[O](22649)'],
    transitionState = 'TS7',
    kinetics = Arrhenius(A=(3293.16,'s^-1'), n=2.5965, Ea=(315.162,'kJ/mol'), T0=(1,'K'), comment="""Estimated using template [R4H_MMS;Cd_rad_out_single;XH_out] for rate rule [R4H_MMS;Cd_rad_out_singleNd;O_H_out]
Euclidian distance = 1.41421356237
family: intra_H_migration"""),
)

reaction(
    label = 'reaction8',
    reactants = ['C=[C][CH]OC#CC=O(25176)'],
    products = ['[CH]=C=COC=C=C[O](22649)'],
    transitionState = 'TS8',
    kinetics = Arrhenius(A=(8.28083e+07,'s^-1'), n=1.36571, Ea=(178.62,'kJ/mol'), T0=(1,'K'), comment="""Estimated using average of templates [RnH;Cd_rad_out_double;Cd_H_out_singleH] + [R5H;Y_rad_out;Cd_H_out_singleH] for rate rule [R5H;Cd_rad_out_double;Cd_H_out_singleH]
Euclidian distance = 2.0
Multiplied by reaction path degeneracy 2.0
family: intra_H_migration"""),
)

reaction(
    label = 'reaction9',
    reactants = ['[CH]=C=[C]OC=C=CO(25177)'],
    products = ['[CH]=C=COC=C=C[O](22649)'],
    transitionState = 'TS9',
    kinetics = Arrhenius(A=(1.45811e+06,'s^-1'), n=1.63609, Ea=(119.541,'kJ/mol'), T0=(1,'K'), comment="""Estimated using average of templates [RnH;Cd_rad_out_double;XH_out] + [R6H;Y_rad_out;XH_out] for rate rule [R6H;Cd_rad_out_double;O_H_out]
Euclidian distance = 2.2360679775
family: intra_H_migration"""),
)

reaction(
    label = 'reaction10',
    reactants = ['[CH]=C=COC=C=C[O](22649)'],
    products = ['C=C=CO[CH][C]=C=O(25178)'],
    transitionState = 'TS10',
    kinetics = Arrhenius(A=(1.86943e+06,'s^-1'), n=1.85754, Ea=(151.521,'kJ/mol'), T0=(1,'K'), comment="""Estimated using template [RnH;Cd_rad_out_singleH;XH_out] for rate rule [R7H;Cd_rad_out_singleH;XH_out]
Euclidian distance = 2.0
family: intra_H_migration"""),
)

reaction(
    label = 'reaction11',
    reactants = ['[CH]=C=C[O](8556)', '[CH]=C=C[O](8556)'],
    products = ['[CH]=C=COC=C=C[O](22649)'],
    transitionState = 'TS11',
    kinetics = Arrhenius(A=(1.63841e+06,'m^3/(mol*s)'), n=0.151, Ea=(0,'kJ/mol'), T0=(1,'K'), comment="""Estimated using template [O_rad/OneDe;Y_rad] for rate rule [O_rad/OneDe;Cd_allenic]
Euclidian distance = 2.0
family: R_Recombination
Ea raised from -0.7 to 0 kJ/mol."""),
)

reaction(
    label = 'reaction12',
    reactants = ['[CH]=C=[CH](18734)', '[O]C=C=C[O](22349)'],
    products = ['[CH]=C=COC=C=C[O](22649)'],
    transitionState = 'TS12',
    kinetics = Arrhenius(A=(6.55363e+06,'m^3/(mol*s)'), n=0.151, Ea=(0,'kJ/mol'), T0=(1,'K'), comment="""Estimated using template [Y_rad;O_rad/OneDe] for rate rule [Cd_allenic;O_rad/OneDe]
Euclidian distance = 2.0
Multiplied by reaction path degeneracy 4.0
family: R_Recombination
Ea raised from -0.7 to 0 kJ/mol."""),
)

reaction(
    label = 'reaction13',
    reactants = ['H(8)', '[CH]=C=CO[C]=C=C[O](25179)'],
    products = ['[CH]=C=COC=C=C[O](22649)'],
    transitionState = 'TS13',
    kinetics = Arrhenius(A=(2.048e+13,'cm^3/(mol*s)'), n=0.206, Ea=(0,'kJ/mol'), T0=(1,'K'), comment="""From training reaction 44 used for Cd_allenic;H_rad
Exact match found for rate rule [Cd_allenic;H_rad]
Euclidian distance = 0
family: R_Recombination
Ea raised from -0.7 to 0 kJ/mol."""),
)

reaction(
    label = 'reaction14',
    reactants = ['H(8)', '[CH]=C=[C]OC=C=C[O](25180)'],
    products = ['[CH]=C=COC=C=C[O](22649)'],
    transitionState = 'TS14',
    kinetics = Arrhenius(A=(2.048e+13,'cm^3/(mol*s)'), n=0.206, Ea=(0,'kJ/mol'), T0=(1,'K'), comment="""From training reaction 44 used for Cd_allenic;H_rad
Exact match found for rate rule [Cd_allenic;H_rad]
Euclidian distance = 0
family: R_Recombination
Ea raised from -0.7 to 0 kJ/mol."""),
)

reaction(
    label = 'reaction15',
    reactants = ['H(8)', '[CH]=C=CO[CH][C]=C=O(25181)'],
    products = ['[CH]=C=COC=C=C[O](22649)'],
    transitionState = 'TS15',
    kinetics = Arrhenius(A=(4.34078e+06,'m^3/(mol*s)'), n=0.278577, Ea=(0,'kJ/mol'), T0=(1,'K'), comment="""Estimated using an average for rate rule [Y_rad;H_rad]
Euclidian distance = 0
family: R_Recombination
Ea raised from -1.4 to 0 kJ/mol."""),
)

reaction(
    label = 'reaction16',
    reactants = ['[CH]=C=COC=C=C[O](22649)'],
    products = ['[CH]=C=COC=C1[CH]O1(25182)'],
    transitionState = 'TS16',
    kinetics = Arrhenius(A=(6.59543e+11,'s^-1'), n=0.400725, Ea=(187.199,'kJ/mol'), T0=(1,'K'), comment="""Estimated using template [R3_D;doublebond_intra;radadd_intra] for rate rule [R3_D;doublebond_intra;radadd_intra_O]
Euclidian distance = 1.0
family: Intra_R_Add_Endocyclic"""),
)

reaction(
    label = 'reaction17',
    reactants = ['[CH]=C=COC=C=C[O](22649)'],
    products = ['[O]C=C=COC1[C]=C1(25183)'],
    transitionState = 'TS17',
    kinetics = Arrhenius(A=(4.2354e+12,'s^-1'), n=-0.1205, Ea=(149.998,'kJ/mol'), T0=(1,'K'), comment="""Estimated using average of templates [Rn;doublebond_intra_CdCdd;radadd_intra_cdsingleH] + [R3;doublebond_intra_CdCdd;radadd_intra_cdsingle] for rate rule [R3;doublebond_intra_CdCdd;radadd_intra_cdsingleH]
Euclidian distance = 1.0
family: Intra_R_Add_Endocyclic"""),
)

reaction(
    label = 'reaction18',
    reactants = ['[CH]=C=COC=C=C[O](22649)'],
    products = ['[CH]=C=COC1[C]=CO1(25184)'],
    transitionState = 'TS18',
    kinetics = Arrhenius(A=(2.47276e+11,'s^-1'), n=0.183155, Ea=(125.658,'kJ/mol'), T0=(1,'K'), comment="""Estimated using average of templates [Rn;doublebond_intra_CdCdd;radadd_intra] + [R4;doublebond_intra;radadd_intra] for rate rule [R4;doublebond_intra_CdCdd;radadd_intra_O]
Euclidian distance = 1.41421356237
family: Intra_R_Add_Endocyclic"""),
)

reaction(
    label = 'reaction19',
    reactants = ['[CH]=C=COC=C=C[O](22649)'],
    products = ['[O]C=C1[CH]OC=C=C1(25114)'],
    transitionState = 'TS19',
    kinetics = Arrhenius(A=(2.16959e+10,'s^-1'), n=0.31, Ea=(12.393,'kJ/mol'), T0=(1,'K'), comment="""Estimated using an average for rate rule [R6_linear;doublebond_intra;radadd_intra_cdsingleH]
Euclidian distance = 0
family: Intra_R_Add_Endocyclic"""),
)

reaction(
    label = 'reaction20',
    reactants = ['[CH]=C=COC=C=C[O](22649)'],
    products = ['[O]C1[C]=COC=C=C1(25185)'],
    transitionState = 'TS20',
    kinetics = Arrhenius(A=(3.22e+12,'s^-1'), n=-0.622, Ea=(148.944,'kJ/mol'), T0=(1,'K'), comment="""From training reaction 832 used for R7;doublebond_intra_CdCdd;radadd_intra_cdsingleH
Exact match found for rate rule [R7;doublebond_intra_CdCdd;radadd_intra_cdsingleH]
Euclidian distance = 0
family: Intra_R_Add_Endocyclic
Ea raised from 144.8 to 148.9 kJ/mol to match endothermicity of reaction."""),
)

reaction(
    label = 'reaction21',
    reactants = ['[CH]=C=COC=C=C[O](22649)'],
    products = ['[CH]=C1[CH]OC=C=CO1(25157)'],
    transitionState = 'TS21',
    kinetics = Arrhenius(A=(8.62196e+06,'s^-1'), n=0.867572, Ea=(62.1704,'kJ/mol'), T0=(1,'K'), comment="""Estimated using template [R7_linear;doublebond_intra;radadd_intra] for rate rule [R7_linear;doublebond_intra;radadd_intra_O]
Euclidian distance = 1.0
family: Intra_R_Add_Endocyclic"""),
)

reaction(
    label = 'reaction22',
    reactants = ['[CH]=C=COC=C=C[O](22649)'],
    products = ['C=C=COC=C=C=O(25186)'],
    transitionState = 'TS22',
    kinetics = Arrhenius(A=(2.1261e+09,'s^-1'), n=0.137, Ea=(24.9733,'kJ/mol'), T0=(1,'K'), comment="""Estimated using template [R7;Y_rad;XH_Rrad] for rate rule [R7radExo;Y_rad;XH_Rrad]
Euclidian distance = 1.0
family: Intra_Disproportionation"""),
)

reaction(
    label = 'reaction23',
    reactants = ['[CH]=C=CO[C]C=C[O](25187)'],
    products = ['[CH]=C=COC=C=C[O](22649)'],
    transitionState = 'TS23',
    kinetics = Arrhenius(A=(2.24409e+10,'s^-1'), n=0.34095, Ea=(22.3009,'kJ/mol'), T0=(1,'K'), comment="""Estimated using average of templates [Rn;Y_rad_NDe;XH_Rrad] + [R2radExo;Y_rad;XH_Rrad] for rate rule [R2radExo;Y_rad_NDe;XH_Rrad]
Euclidian distance = 1.0
family: Intra_Disproportionation"""),
)

reaction(
    label = 'reaction24',
    reactants = ['[CH]=C=CO[CH]C=[C][O](25188)'],
    products = ['[CH]=C=COC=C=C[O](22649)'],
    transitionState = 'TS24',
    kinetics = Arrhenius(A=(1.949e+11,'s^-1'), n=0.486, Ea=(22.8614,'kJ/mol'), T0=(1,'K'), comment="""Estimated using an average for rate rule [R2radExo;Y_rad;XH_Rrad_NDe]
Euclidian distance = 0
family: Intra_Disproportionation"""),
)

reaction(
    label = 'reaction25',
    reactants = ['[CH]=[C]CO[C]=C=C[O](25189)'],
    products = ['[CH]=C=COC=C=C[O](22649)'],
    transitionState = 'TS25',
    kinetics = Arrhenius(A=(5.2748e+09,'s^-1'), n=0.37, Ea=(88.9686,'kJ/mol'), T0=(1,'K'), comment="""Estimated using average of templates [R3;Y_rad_De;XH_Rrad] + [R3radExo;Y_rad;XH_Rrad] for rate rule [R3radExo;Y_rad_De;XH_Rrad]
Euclidian distance = 1.0
Multiplied by reaction path degeneracy 2.0
family: Intra_Disproportionation"""),
)

reaction(
    label = 'reaction26',
    reactants = ['[CH]=C=[C]OC[C]=C[O](25190)'],
    products = ['[CH]=C=COC=C=C[O](22649)'],
    transitionState = 'TS26',
    kinetics = Arrhenius(A=(1.08e+10,'s^-1'), n=-0.305, Ea=(93.094,'kJ/mol'), T0=(1,'K'), Tmin=(300,'K'), Tmax=(1500,'K'), comment="""Estimated using template [R3;Y_rad_De;XH_Rrad_De] for rate rule [R3radExo;Y_rad_De;XH_Rrad_De]
Euclidian distance = 1.0
Multiplied by reaction path degeneracy 2.0
family: Intra_Disproportionation"""),
)

reaction(
    label = 'reaction27',
    reactants = ['[CH]=C[CH]O[C]=C=C[O](25191)'],
    products = ['[CH]=C=COC=C=C[O](22649)'],
    transitionState = 'TS27',
    kinetics = Arrhenius(A=(4.48312e+08,'s^-1'), n=0.311, Ea=(39.225,'kJ/mol'), T0=(1,'K'), comment="""Estimated using template [R4;Y_rad_De;XH_Rrad] for rate rule [R4radEndo;Y_rad_De;XH_Rrad]
Euclidian distance = 1.0
family: Intra_Disproportionation"""),
)

reaction(
    label = 'reaction28',
    reactants = ['[CH]=C=[C]O[CH]C=C[O](25192)'],
    products = ['[CH]=C=COC=C=C[O](22649)'],
    transitionState = 'TS28',
    kinetics = Arrhenius(A=(4.48312e+08,'s^-1'), n=0.311, Ea=(39.225,'kJ/mol'), T0=(1,'K'), comment="""Estimated using template [R4;Y_rad_De;XH_Rrad] for rate rule [R4radEndo;Y_rad_De;XH_Rrad]
Euclidian distance = 1.0
family: Intra_Disproportionation"""),
)

reaction(
    label = 'reaction29',
    reactants = ['[CH]=C=[C]OC=[C]C[O](25193)'],
    products = ['[CH]=C=COC=C=C[O](22649)'],
    transitionState = 'TS29',
    kinetics = Arrhenius(A=(3.70659e+09,'s^-1'), n=0.137, Ea=(39.225,'kJ/mol'), T0=(1,'K'), comment="""Estimated using template [R5;Y_rad_De;XH_Rrad] for rate rule [R5radEndo;Y_rad_De;XH_Rrad]
Euclidian distance = 1.0
Multiplied by reaction path degeneracy 2.0
family: Intra_Disproportionation"""),
)

reaction(
    label = 'reaction30',
    reactants = ['[CH]=[C]CO[CH][C]=C=O(25194)'],
    products = ['[CH]=C=COC=C=C[O](22649)'],
    transitionState = 'TS30',
    kinetics = Arrhenius(A=(4.25221e+09,'s^-1'), n=0.137, Ea=(24.9733,'kJ/mol'), T0=(1,'K'), comment="""Estimated using template [R5;Y_rad;XH_Rrad] for rate rule [R5radExo;Y_rad;XH_Rrad]
Euclidian distance = 1.0
Multiplied by reaction path degeneracy 2.0
family: Intra_Disproportionation"""),
)

reaction(
    label = 'reaction31',
    reactants = ['[CH]C=CO[CH][C]=C=O(25195)'],
    products = ['[CH]=C=COC=C=C[O](22649)'],
    transitionState = 'TS31',
    kinetics = Arrhenius(A=(2.1261e+09,'s^-1'), n=0.137, Ea=(24.9733,'kJ/mol'), T0=(1,'K'), comment="""Estimated using template [R6;Y_rad;XH_Rrad] for rate rule [R6radEndo;Y_rad;XH_Rrad]
Euclidian distance = 1.0
family: Intra_Disproportionation"""),
)

reaction(
    label = 'reaction32',
    reactants = ['[CH]=C=COC=C=C[O](22649)'],
    products = ['[CH]=C=C[C](C=O)C=O(25196)'],
    transitionState = 'TS32',
    kinetics = Arrhenius(A=(7040,'s^-1'), n=2.66, Ea=(313.8,'kJ/mol'), T0=(1,'K'), Tmin=(300,'K'), Tmax=(1500,'K'), comment="""Estimated using an average for rate rule [R_ROR;R1_doublebond;R2_doublebond_H;R_O_C]
Euclidian distance = 0
family: ketoenol"""),
)

reaction(
    label = 'reaction36',
    reactants = ['[CH]=C=COC=C=C[O](22649)'],
    products = ['[CH]=C(C=O)C=C=C[O](22652)'],
    transitionState = 'TS33',
    kinetics = Arrhenius(A=(7040,'s^-1'), n=2.66, Ea=(313.8,'kJ/mol'), T0=(1,'K'), Tmin=(300,'K'), Tmax=(1500,'K'), comment="""Estimated using an average for rate rule [R_ROR;R1_doublebond;R2_doublebond_H;R_O_C]
Euclidian distance = 0
family: ketoenol"""),
)

reaction(
    label = 'reaction34',
    reactants = ['[CH]C1=CO[CH]C1=C[O](25197)'],
    products = ['[CH]=C=COC=C=C[O](22649)'],
    transitionState = 'TS34',
    kinetics = Arrhenius(A=(7.69248e+14,'s^-1'), n=-0.917475, Ea=(150.312,'kJ/mol'), T0=(1,'K'), Tmin=(303.03,'K'), Tmax=(2000,'K'), comment="""From training reaction 1 used for R5JJ
Exact match found for rate rule [R5JJ]
Euclidian distance = 0
family: 1,4_Cyclic_birad_scission"""),
)

reaction(
    label = 'reaction35',
    reactants = ['O(T)(63)', '[CH]=C=COC=C=[CH](22632)'],
    products = ['[CH]=C=COC=C=C[O](22649)'],
    transitionState = 'TS35',
    kinetics = Arrhenius(A=(187219,'m^3/(mol*s)'), n=1.13083, Ea=(163.847,'kJ/mol'), T0=(1,'K'), Tmin=(303.03,'K'), Tmax=(2000,'K'), comment="""Estimated using template [Y_rad;O_birad] for rate rule [Cd_pri_rad;O_birad]
Euclidian distance = 2.0
Multiplied by reaction path degeneracy 2.0
family: Birad_R_Recombination"""),
)

reaction(
    label = 'reaction36',
    reactants = ['H(8)', '[C]#C[CH]OC=C=C[O](25198)'],
    products = ['[CH]=C=COC=C=C[O](22649)'],
    transitionState = 'TS36',
    kinetics = Arrhenius(A=(1e+07,'m^3/(mol*s)'), n=0, Ea=(0,'kJ/mol'), T0=(1,'K'), comment="""Estimated using an average for rate rule [H_rad;Birad]
Euclidian distance = 0
family: Birad_R_Recombination"""),
)

reaction(
    label = 'reaction37',
    reactants = ['[CH]=C=COC=C=C[O](22649)'],
    products = ['C#CC1OC1[C]=C[O](25199)'],
    transitionState = 'TS37',
    kinetics = Arrhenius(A=(1e+10,'s^-1'), n=0, Ea=(108.156,'kJ/mol'), T0=(1,'K'), comment="""Estimated using template [R4_S_D;doublebond_intra;radadd_intra_csHDe] for rate rule [R4_S_D;doublebond_intra;radadd_intra_csHCt]
Euclidian distance = 1.0
family: Intra_R_Add_Exocyclic"""),
)

reaction(
    label = 'reaction38',
    reactants = ['[CH]=C=C[O](8556)', 'C#CC=O(21959)'],
    products = ['[CH]=C=COC=C=C[O](22649)'],
    transitionState = 'TS38',
    kinetics = Arrhenius(A=(1.07e+06,'cm^3/(mol*s)'), n=2.43, Ea=(22.5936,'kJ/mol'), T0=(1,'K'), Tmin=(300,'K'), Tmax=(1500,'K'), comment="""Estimated using template [Od_CO-CtH;YJ] for rate rule [Od_CO-CtH;CdsJ=Cdd]
Euclidian distance = 2.0
family: R_Addition_MultipleBond"""),
)

reaction(
    label = 'reaction39',
    reactants = ['C#CCO[C]=C=C[O](25200)'],
    products = ['[CH]=C=COC=C=C[O](22649)'],
    transitionState = 'TS39',
    kinetics = Arrhenius(A=(3.32e+08,'s^-1'), n=1.34, Ea=(150.415,'kJ/mol'), T0=(1,'K'), comment="""Estimated using average of templates [R3H_SS;Cd_rad_out_double;Cs_H_out_1H] + [R3H_SS_O;Y_rad_out;Cs_H_out_H/OneDe] for rate rule [R3H_SS_O;Cd_rad_out_double;Cs_H_out_H/Ct]
Euclidian distance = 2.2360679775
Multiplied by reaction path degeneracy 2.0
family: intra_H_migration"""),
)

reaction(
    label = 'reaction40',
    reactants = ['[C]#CCOC=C=C[O](25201)'],
    products = ['[CH]=C=COC=C=C[O](22649)'],
    transitionState = 'TS40',
    kinetics = Arrhenius(A=(1.47715e+10,'s^-1'), n=0.8, Ea=(147.277,'kJ/mol'), T0=(1,'K'), comment="""Estimated using template [R3H_MS;Y_rad_out;Cs_H_out_1H] for rate rule [R3H_TS;Ct_rad_out;Cs_H_out_H/NonDeO]
Euclidian distance = 1.73205080757
Multiplied by reaction path degeneracy 2.0
family: intra_H_migration"""),
)

reaction(
    label = 'reaction41',
    reactants = ['C#CCO[CH][C]=C=O(25202)'],
    products = ['[CH]=C=COC=C=C[O](22649)'],
    transitionState = 'TS41',
    kinetics = Arrhenius(A=(59826.2,'s^-1'), n=1.86494, Ea=(62.4776,'kJ/mol'), T0=(1,'K'), comment="""Estimated using template [R5H;Y_rad_out;Cs_H_out_H/OneDe] for rate rule [R5H;Cd_rad_out;Cs_H_out_H/Ct]
Euclidian distance = 1.41421356237
Multiplied by reaction path degeneracy 2.0
family: intra_H_migration"""),
)

reaction(
    label = 'reaction42',
    reactants = ['[C]#C[CH]OC=C=CO(25203)'],
    products = ['[CH]=C=COC=C=C[O](22649)'],
    transitionState = 'TS42',
    kinetics = Arrhenius(A=(3728.31,'s^-1'), n=2.31462, Ea=(124.436,'kJ/mol'), T0=(1,'K'), comment="""Estimated using average of templates [RnH;Y_rad_out;O_H_out] + [R8Hall;Y_rad_out;XH_out] for rate rule [R8Hall;Ct_rad_out;O_H_out]
Euclidian distance = 1.41421356237
family: intra_H_migration"""),
)

reaction(
    label = 'reaction43',
    reactants = ['[CH]=C=COC=C=C[O](22649)'],
    products = ['C#CC1O[CH]C1=C[O](25204)'],
    transitionState = 'TS43',
    kinetics = Arrhenius(A=(8.33681e+07,'s^-1'), n=1.1093, Ea=(126.639,'kJ/mol'), T0=(1,'K'), comment="""Estimated using template [R4_S_D;doublebond_intra;radadd_intra_cs] for rate rule [R4_S_D;doublebond_intra;radadd_intra_csHCt]
Euclidian distance = 2.0
family: Intra_R_Add_Endocyclic"""),
)

reaction(
    label = 'reaction44',
    reactants = ['[CH]=C=COC=C=C[O](22649)'],
    products = ['C#CC1OC=[C]C1[O](25205)'],
    transitionState = 'TS44',
    kinetics = Arrhenius(A=(8.68775e+12,'s^-1'), n=-0.482695, Ea=(142.968,'kJ/mol'), T0=(1,'K'), comment="""Estimated using template [R5_linear;doublebond_intra;radadd_intra_csHDe] for rate rule [R5_linear;doublebond_intra_CdCdd;radadd_intra_csHCt]
Euclidian distance = 1.41421356237
family: Intra_R_Add_Endocyclic"""),
)

reaction(
    label = 'reaction45',
    reactants = ['[CH]=C=COC=C=C[O](22649)'],
    products = ['[C]1[CH]OC=C=COC=1(25206)'],
    transitionState = 'TS45',
    kinetics = Arrhenius(A=(8.28217e+09,'s^-1'), n=0.356567, Ea=(65.4205,'kJ/mol'), T0=(1,'K'), comment="""Estimated using template [R6plus;triplebond_intra_H;radadd_intra] for rate rule [R8_linear;triplebond_intra_H;radadd_intra_O]
Euclidian distance = 2.2360679775
family: Intra_R_Add_Endocyclic"""),
)

reaction(
    label = 'reaction46',
    reactants = ['[CH]=C=COC=C=C[O](22649)'],
    products = ['C#CCOC=C=C=O(25207)'],
    transitionState = 'TS46',
    kinetics = Arrhenius(A=(1.85329e+09,'s^-1'), n=0.137, Ea=(39.225,'kJ/mol'), T0=(1,'K'), comment="""Estimated using template [R5;Y_rad_De;XH_Rrad] for rate rule [R5radExo;Y_rad_De;XH_Rrad]
Euclidian distance = 1.0
family: Intra_Disproportionation"""),
)

reaction(
    label = 'reaction47',
    reactants = ['[C]#C(5143)', '[CH]OC=C=C[O](23145)'],
    products = ['[CH]=C=COC=C=C[O](22649)'],
    transitionState = 'TS47',
    kinetics = Arrhenius(A=(1.14854e+06,'m^3/(mol*s)'), n=0.575199, Ea=(34.3157,'kJ/mol'), T0=(1,'K'), comment="""Estimated using template [Y_rad;Birad] for rate rule [Ct_rad;Birad]
Euclidian distance = 1.0
family: Birad_R_Recombination"""),
)

reaction(
    label = 'reaction48',
    reactants = ['H(8)', '[CH]=C=COC#CC=O(25208)'],
    products = ['[CH]=C=COC=C=C[O](22649)'],
    transitionState = 'TS48',
    kinetics = Arrhenius(A=(4278.27,'m^3/(mol*s)'), n=1.383, Ea=(15.1097,'kJ/mol'), T0=(1,'K'), comment="""Estimated using template [Ct_Ct;HJ] for rate rule [Ct-O_Ct;HJ]
Euclidian distance = 1.0
family: R_Addition_MultipleBond"""),
)

reaction(
    label = 'reaction49',
    reactants = ['[CH]=C=CO[C]=CC=O(25209)'],
    products = ['[CH]=C=COC=C=C[O](22649)'],
    transitionState = 'TS49',
    kinetics = Arrhenius(A=(1.75e+11,'s^-1'), n=0.633, Ea=(196.46,'kJ/mol'), T0=(1,'K'), comment="""Estimated using an average for rate rule [R2H_D;Cd_rad_out_singleNd;Cd_H_out_singleDe]
Euclidian distance = 0
family: intra_H_migration"""),
)

reaction(
    label = 'reaction50',
    reactants = ['[CH]=C=COC=C=C[O](22649)'],
    products = ['[CH]=C=CO[CH]C=C=O(25210)'],
    transitionState = 'TS50',
    kinetics = Arrhenius(A=(4.96519e+09,'s^-1'), n=1.05826, Ea=(162.779,'kJ/mol'), T0=(1,'K'), comment="""Estimated using template [R2H_S;Cd_rad_out_Cd;XH_out] for rate rule [R2H_S;Cd_rad_out_Cd;CO_H_out]
Euclidian distance = 1.0
family: intra_H_migration"""),
)

reaction(
    label = 'reaction51',
    reactants = ['[CH]=C=[C]OC=CC=O(25211)'],
    products = ['[CH]=C=COC=C=C[O](22649)'],
    transitionState = 'TS51',
    kinetics = Arrhenius(A=(9.01194e+11,'s^-1'), n=1.09397, Ea=(394.305,'kJ/mol'), T0=(1,'K'), comment="""Estimated using template [R4H_RSD;Cd_rad_out;Cd_H_out_single] for rate rule [R4H_SSD;Cd_rad_out_double;Cd_H_out_singleDe]
Euclidian distance = 1.73205080757
family: intra_H_migration"""),
)

reaction(
    label = 'reaction52',
    reactants = ['[CH]=C=COC=C=C[O](22649)'],
    products = ['[CH]C1=COC=C1C=O(25071)'],
    transitionState = 'TS52',
    kinetics = Arrhenius(A=(2.42998e+09,'s^-1'), n=0.57, Ea=(38.3882,'kJ/mol'), T0=(1,'K'), comment="""Estimated using template [R5_DS_D;doublebond_intra;radadd_intra_cdsingle] for rate rule [R5_DS_D;doublebond_intra;radadd_intra_cdsingleDe]
Euclidian distance = 1.0
family: Intra_R_Add_Endocyclic"""),
)

reaction(
    label = 'reaction53',
    reactants = ['[CH]=O(373)', '[C]=COC=C=[CH](23405)'],
    products = ['[CH]=C=COC=C=C[O](22649)'],
    transitionState = 'TS53',
    kinetics = Arrhenius(A=(1.14854e+06,'m^3/(mol*s)'), n=0.575199, Ea=(34.3157,'kJ/mol'), T0=(1,'K'), comment="""Estimated using template [Y_rad;Birad] for rate rule [CO_pri_rad;Birad]
Euclidian distance = 2.0
family: Birad_R_Recombination"""),
)

reaction(
    label = 'reaction54',
    reactants = ['[C]#C[CH]OC=CC=O(25212)'],
    products = ['[CH]=C=COC=C=C[O](22649)'],
    transitionState = 'TS54',
    kinetics = Arrhenius(A=(2.67387e+08,'s^-1'), n=1.29644, Ea=(188.354,'kJ/mol'), T0=(1,'K'), comment="""Estimated using template [RnH;Y_rad_out;Cd_H_out_singleDe] for rate rule [R6HJ_2;Ct_rad_out;Cd_H_out_singleDe]
Euclidian distance = 2.2360679775
family: intra_H_migration"""),
)

reaction(
    label = 'reaction55',
    reactants = ['[CH]=C=COC=C=C[O](22649)'],
    products = ['C#CC1OC=C1C=O(25213)'],
    transitionState = 'TS55',
    kinetics = Arrhenius(A=(2e+12,'s^-1'), n=0, Ea=(7.5312,'kJ/mol'), T0=(1,'K'), comment="""Estimated using template [Rn;C_rad_out_H/OneDe;Ypri_rad_out] for rate rule [R4_SSD;C_rad_out_H/OneDe;CdsinglepriDe_rad_out]
Euclidian distance = 2.82842712475
family: Birad_recombination"""),
)

network(
    label = '4952',
    isomers = [
        '[CH]=C=COC=C=C[O](22649)',
    ],
    reactants = [
        ('C#CC=O(21959)', 'C#CC=O(21959)'),
    ],
    bathGas = {
        'N2': 0.25,
        'Ne': 0.25,
        'He': 0.25,
        'Ar': 0.25,
    },
)

pressureDependence(
    label = '4952',
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

