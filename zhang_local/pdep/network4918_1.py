species(
    label = '[CH]C(C=O)=C[C]=C(22615)',
    structure = SMILES('[CH]C(C=O)=C[C]=C'),
    E0 = (496.77,'kJ/mol'),
    modes = [
        HarmonicOscillator(frequencies=([1685,370,2782.5,750,1395,475,1775,1000,2950,3100,1380,975,1025,1650,350,440,435,1725,3010,987.5,1337.5,450,1655,180,180,180,180],'cm^-1')),
        HinderedRotor(inertia=(2.09979,'amu*angstrom^2'), symmetry=1, barrier=(48.2784,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(2.11073,'amu*angstrom^2'), symmetry=1, barrier=(48.5298,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(2.1197,'amu*angstrom^2'), symmetry=1, barrier=(48.7362,'kJ/mol'), semiclassical=False),
    ],
    spinMultiplicity = 4,
    opticalIsomers = 1,
    molecularWeight = (93.1033,'amu'),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[1.29383,0.062448,-5.49929e-05,2.6954e-08,-5.66625e-12,59842.6,22.8521], Tmin=(100,'K'), Tmax=(1098.96,'K')), NASAPolynomial(coeffs=[9.25525,0.0334701,-1.54402e-05,2.96006e-09,-2.07929e-13,58092.7,-16.3084], Tmin=(1098.96,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(496.77,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(270.22,'J/(mol*K)'), comment="""Thermo group additivity estimation: group(Cs-(Cds-Cds)HHH) + group(Cd-CdCs(CO)) + group(Cds-Cds(Cds-Cds)H) + group(Cds-Cds(Cds-Cds)H) + group(Cds-O2d(Cds-Cds)H) + group(Cds-CdsHH) + radical(C=CJC=C) + radical(AllylJ2_triplet)"""),
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
    label = 'C#C[CH2](17441)',
    structure = SMILES('C#C[CH2]'),
    E0 = (328.481,'kJ/mol'),
    modes = [
        HarmonicOscillator(frequencies=([3000,3100,440,815,1455,1000,2175,525,1131.03,1132.16,1135.9],'cm^-1')),
        HinderedRotor(inertia=(0.154206,'amu*angstrom^2'), symmetry=1, barrier=(3.5455,'kJ/mol'), semiclassical=False),
    ],
    spinMultiplicity = 2,
    opticalIsomers = 1,
    molecularWeight = (39.0559,'amu'),
    collisionModel = TransportData(shapeIndex=2, epsilon=(2095.25,'J/mol'), sigma=(4.76,'angstroms'), dipoleMoment=(0,'C*m'), polarizability=(0,'angstroms^3'), rotrelaxcollnum=1.0, comment="""GRI-Mech"""),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[3.32026,0.0108736,8.62061e-06,-1.82973e-08,7.68649e-12,39535.3,8.27851], Tmin=(100,'K'), Tmax=(960.555,'K')), NASAPolynomial(coeffs=[6.38511,0.00814486,-2.78734e-06,4.95348e-10,-3.50148e-14,38483.6,-8.79383], Tmin=(960.555,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(328.481,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(128.874,'J/(mol*K)'), comment="""Thermo group additivity estimation: group(Cs-CtHHH) + group(Ct-CtCs) + group(Ct-CtH) + radical(Propargyl)"""),
)

species(
    label = '[CH]=C1[CH]C(=C)C1[O](24927)',
    structure = SMILES('[CH]=C1[CH]C(=C)C1[O]'),
    E0 = (597.485,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
    molecularWeight = (93.1033,'amu'),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[2.05283,0.0184118,8.32106e-05,-1.29788e-07,5.32125e-11,71953.3,22.9707], Tmin=(100,'K'), Tmax=(940.36,'K')), NASAPolynomial(coeffs=[18.8488,0.00909785,-1.03951e-06,2.03187e-10,-2.57447e-14,66047.4,-71.6337], Tmin=(940.36,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(597.485,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(282.692,'J/(mol*K)'), comment="""Thermo group additivity estimation: group(O2s-CsH) + group(Cs-(Cds-Cds)(Cds-Cds)OsH) + group(Cs-(Cds-Cds)(Cds-Cds)HH) + group(Cds-CdsCsCs) + group(Cds-CdsCsCs) + group(Cds-CdsHH) + group(Cds-CdsHH) + ring(Cyclobutane) + radical(Cds_P) + radical(C=CCJC=C) + radical(CC(C)OJ)"""),
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
    label = '[CH]=C(C#C[CH2])C=O(24935)',
    structure = SMILES('[CH]=C(C#C[CH2])C=O'),
    E0 = (518.076,'kJ/mol'),
    modes = [
        HarmonicOscillator(frequencies=([3120,650,792.5,1650,2782.5,750,1395,475,1775,1000,350,440,435,1725,2100,2250,500,550,3000,3100,440,815,1455,1000],'cm^-1')),
        HinderedRotor(inertia=(1.24139,'amu*angstrom^2'), symmetry=1, barrier=(28.5421,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(1.24179,'amu*angstrom^2'), symmetry=1, barrier=(28.5513,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(1.2417,'amu*angstrom^2'), symmetry=1, barrier=(28.549,'kJ/mol'), semiclassical=False),
    ],
    spinMultiplicity = 3,
    opticalIsomers = 1,
    molecularWeight = (92.0954,'amu'),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[1.08045,0.0571037,-5.13689e-05,2.19254e-08,-3.65789e-12,62421.1,23.4975], Tmin=(100,'K'), Tmax=(1446.05,'K')), NASAPolynomial(coeffs=[16.6694,0.0139823,-6.63877e-06,1.30367e-09,-9.27002e-14,57912.6,-57.4598], Tmin=(1446.05,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(518.076,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(245.277,'J/(mol*K)'), comment="""Thermo group additivity estimation: group(Cs-CtHHH) + group(Cds-Cds(Cds-O2d)(Cds-Cds)) + group(Cds-O2d(Cds-Cds)H) + group(Cds-CdsHH) + group(Ct-CtCs) + group(Ct-Ct(Cds-Cds)) + radical(Cds_P) + radical(Propargyl)"""),
)

species(
    label = '[CH]C(C=O)=CC#C(22634)',
    structure = SMILES('[CH]C(C=O)=CC#C'),
    E0 = (474.56,'kJ/mol'),
    modes = [
        HarmonicOscillator(frequencies=([2782.5,750,1395,475,1775,1000,750,770,3400,2100,350,440,435,1725,2175,525,3010,987.5,1337.5,450,1655,180,180,180],'cm^-1')),
        HinderedRotor(inertia=(2.0549,'amu*angstrom^2'), symmetry=1, barrier=(47.2462,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(2.05334,'amu*angstrom^2'), symmetry=1, barrier=(47.2103,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(2.0545,'amu*angstrom^2'), symmetry=1, barrier=(47.237,'kJ/mol'), semiclassical=False),
    ],
    spinMultiplicity = 3,
    opticalIsomers = 1,
    molecularWeight = (92.0954,'amu'),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[1.52185,0.0555926,-4.20631e-05,1.56625e-08,-2.41268e-12,57164.2,21.3186], Tmin=(100,'K'), Tmax=(1465.47,'K')), NASAPolynomial(coeffs=[11.8146,0.0274985,-1.33072e-05,2.58108e-09,-1.81074e-13,54147.4,-32.272], Tmin=(1465.47,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(474.56,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(245.277,'J/(mol*K)'), comment="""Thermo group additivity estimation: group(Cs-(Cds-Cds)HHH) + group(Cd-CdCs(CO)) + group(Cds-CdsCtH) + group(Cds-O2d(Cds-Cds)H) + group(Ct-Ct(Cds-Cds)) + group(Ct-CtH) + radical(AllylJ2_triplet)"""),
)

species(
    label = '[CH]C(=[C]C=C)C=O(24936)',
    structure = SMILES('[CH]C(=[C]C=C)C=O'),
    E0 = (496.77,'kJ/mol'),
    modes = [
        HarmonicOscillator(frequencies=([1685,370,2782.5,750,1395,475,1775,1000,2950,3100,1380,975,1025,1650,350,440,435,1725,3010,987.5,1337.5,450,1655,180,180,180,180],'cm^-1')),
        HinderedRotor(inertia=(2.09979,'amu*angstrom^2'), symmetry=1, barrier=(48.2784,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(2.11073,'amu*angstrom^2'), symmetry=1, barrier=(48.5298,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(2.1197,'amu*angstrom^2'), symmetry=1, barrier=(48.7362,'kJ/mol'), semiclassical=False),
    ],
    spinMultiplicity = 4,
    opticalIsomers = 1,
    molecularWeight = (93.1033,'amu'),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[1.29383,0.062448,-5.49929e-05,2.6954e-08,-5.66625e-12,59842.6,22.8521], Tmin=(100,'K'), Tmax=(1098.96,'K')), NASAPolynomial(coeffs=[9.25525,0.0334701,-1.54402e-05,2.96006e-09,-2.07929e-13,58092.7,-16.3084], Tmin=(1098.96,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(496.77,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(270.22,'J/(mol*K)'), comment="""Thermo group additivity estimation: group(Cs-(Cds-Cds)HHH) + group(Cd-CdCs(CO)) + group(Cds-Cds(Cds-Cds)H) + group(Cds-Cds(Cds-Cds)H) + group(Cds-O2d(Cds-Cds)H) + group(Cds-CdsHH) + radical(AllylJ2_triplet) + radical(C=CJC=C)"""),
)

species(
    label = '[CH]C(C=O)=CC=[CH](24771)',
    structure = SMILES('[CH]C(C=O)=CC=[CH]'),
    E0 = (544.87,'kJ/mol'),
    modes = [
        HarmonicOscillator(frequencies=([3120,650,792.5,1650,2782.5,750,1395,475,1775,1000,350,440,435,1725,2995,3025,975,1000,1300,1375,400,500,1630,1680,180,180,180],'cm^-1')),
        HinderedRotor(inertia=(2.04874,'amu*angstrom^2'), symmetry=1, barrier=(47.1045,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(2.04903,'amu*angstrom^2'), symmetry=1, barrier=(47.1112,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(2.04866,'amu*angstrom^2'), symmetry=1, barrier=(47.1027,'kJ/mol'), semiclassical=False),
    ],
    spinMultiplicity = 4,
    opticalIsomers = 1,
    molecularWeight = (93.1033,'amu'),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[1.15442,0.061202,-4.89738e-05,1.98015e-08,-3.28945e-12,65636,24.151], Tmin=(100,'K'), Tmax=(1392.79,'K')), NASAPolynomial(coeffs=[12.8983,0.027474,-1.26491e-05,2.41419e-09,-1.68458e-13,62364.7,-36.397], Tmin=(1392.79,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(544.87,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(270.22,'J/(mol*K)'), comment="""Thermo group additivity estimation: group(Cs-(Cds-Cds)HHH) + group(Cd-CdCs(CO)) + group(Cds-Cds(Cds-Cds)H) + group(Cds-Cds(Cds-Cds)H) + group(Cds-O2d(Cds-Cds)H) + group(Cds-CdsHH) + radical(AllylJ2_triplet) + radical(Cds_P)"""),
)

species(
    label = '[CH]C(=C=O)[CH]C=C(24937)',
    structure = SMILES('[CH]C(=C=O)[CH]C=C'),
    E0 = (443.645,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
    molecularWeight = (93.1033,'amu'),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[1.19625,0.0587673,-4.55133e-05,1.87068e-08,-3.19842e-12,53461.4,22.8321], Tmin=(100,'K'), Tmax=(1355.61,'K')), NASAPolynomial(coeffs=[11.5017,0.028359,-1.18661e-05,2.15965e-09,-1.46806e-13,50667.4,-30.021], Tmin=(1355.61,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(443.645,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(270.22,'J/(mol*K)'), comment="""Thermo group additivity estimation: group(Cs-Cd(CCO)HH) + group(Cs-(Cds-Cdd-O2d)HHH) + group(Cds-(Cdd-O2d)CsCs) + other(ketene_1C-C_1C-H) + group(Cds-CdsCsH) + group(Cds-CdsHH) + radical(AllylJ2_triplet) + radical(C=CCJC=C=O)"""),
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
    label = '[CH]C(C#C[CH2])=C[O](24938)',
    structure = SMILES('[CH]C(C#C[CH2])=C[O]'),
    E0 = (647.832,'kJ/mol'),
    modes = [
        HarmonicOscillator(frequencies=([3000,3100,440,815,1455,1000,350,440,435,1725,2100,2250,500,550,3010,987.5,1337.5,450,1655,314.431,314.431,314.432,314.433,314.436],'cm^-1')),
        HinderedRotor(inertia=(0.707442,'amu*angstrom^2'), symmetry=1, barrier=(49.6326,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.707432,'amu*angstrom^2'), symmetry=1, barrier=(49.6327,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.707423,'amu*angstrom^2'), symmetry=1, barrier=(49.6327,'kJ/mol'), semiclassical=False),
    ],
    spinMultiplicity = 5,
    opticalIsomers = 1,
    molecularWeight = (92.0954,'amu'),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[1.1702,0.0515381,-1.88326e-05,-1.78464e-08,1.21153e-11,78027.7,25.4748], Tmin=(100,'K'), Tmax=(957.505,'K')), NASAPolynomial(coeffs=[15.3028,0.0179605,-6.1184e-06,1.07321e-09,-7.55795e-14,74154.1,-48.1882], Tmin=(957.505,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(647.832,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(245.277,'J/(mol*K)'), comment="""Thermo group additivity estimation: group(O2s-(Cds-Cd)H) + group(Cs-(Cds-Cds)HHH) + group(Cs-CtHHH) + group(Cds-CdsCtCs) + group(Cds-CdsOsH) + group(Ct-CtCs) + group(Ct-Ct(Cds-Cds)) + radical(AllylJ2_triplet) + radical(Propargyl) + radical(C=COJ)"""),
)

species(
    label = '[CH]C([C]=O)=C[C]=C(24939)',
    structure = SMILES('[CH]C([C]=O)=C[C]=C'),
    E0 = (660.373,'kJ/mol'),
    modes = [
        HarmonicOscillator(frequencies=([1685,370,2950,3100,1380,975,1025,1650,350,440,435,1725,1855,455,950,3010,987.5,1337.5,450,1655,180,180,180,180],'cm^-1')),
        HinderedRotor(inertia=(2.09695,'amu*angstrom^2'), symmetry=1, barrier=(48.2131,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(2.09751,'amu*angstrom^2'), symmetry=1, barrier=(48.2258,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(2.09621,'amu*angstrom^2'), symmetry=1, barrier=(48.196,'kJ/mol'), semiclassical=False),
    ],
    spinMultiplicity = 5,
    opticalIsomers = 1,
    molecularWeight = (92.0954,'amu'),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[1.17766,0.0662367,-7.31014e-05,4.61521e-08,-1.22628e-11,79522.7,23.8384], Tmin=(100,'K'), Tmax=(897.248,'K')), NASAPolynomial(coeffs=[8.99464,0.031388,-1.48421e-05,2.86465e-09,-2.01685e-13,78119.9,-13.0265], Tmin=(897.248,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(660.373,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(245.277,'J/(mol*K)'), comment="""Thermo group additivity estimation: group(Cs-(Cds-Cds)HHH) + group(Cd-CdCs(CO)) + group(Cds-Cds(Cds-Cds)H) + group(Cds-Cds(Cds-Cds)H) + group(Cds-O2d(Cds-Cds)H) + group(Cds-CdsHH) + radical(AllylJ2_triplet) + radical(C=CJC=C) + radical(C=C(C)CJ=O)"""),
)

species(
    label = '[CH]C(C=C=[CH])=C[O](22639)',
    structure = SMILES('[CH]C(C=C=[CH])=C[O]'),
    E0 = (654.74,'kJ/mol'),
    modes = [
        HarmonicOscillator(frequencies=([3120,650,792.5,1650,540,610,2055,350,440,435,1725,2995,3025,975,1000,1300,1375,400,500,1630,1680,180,180,180,180],'cm^-1')),
        HinderedRotor(inertia=(2.0303,'amu*angstrom^2'), symmetry=1, barrier=(46.6806,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(2.03052,'amu*angstrom^2'), symmetry=1, barrier=(46.6856,'kJ/mol'), semiclassical=False),
    ],
    spinMultiplicity = 5,
    opticalIsomers = 1,
    molecularWeight = (92.0954,'amu'),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[-0.459934,0.0743252,-7.39542e-05,3.5466e-08,-6.30027e-12,78928.1,28.4852], Tmin=(100,'K'), Tmax=(1631.97,'K')), NASAPolynomial(coeffs=[19.9463,0.009355,-4.93136e-07,-1.57955e-10,1.70616e-14,74259,-73.857], Tmin=(1631.97,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(654.74,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(249.434,'J/(mol*K)'), comment="""Thermo group additivity estimation: group(O2s-(Cds-Cd)H) + group(Cs-(Cds-Cds)HHH) + group(Cds-Cds(Cds-Cds)Cs) + group(Cds-Cds(Cds-Cds)H) + group(Cds-CdsOsH) + group(Cds-CdsHH) + group(Cdd-CdsCds) + radical(C=C=CJ) + radical(AllylJ2_triplet) + radical(C=COJ)"""),
)

species(
    label = '[CH]C1([CH]C1=C)C=O(24940)',
    structure = SMILES('[CH]C1([CH]C1=C)C=O'),
    E0 = (685.636,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
    molecularWeight = (93.1033,'amu'),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[1.30506,0.0501309,-2.63387e-05,-5.59805e-09,6.73976e-12,82568.1,24.599], Tmin=(100,'K'), Tmax=(990.714,'K')), NASAPolynomial(coeffs=[14.4765,0.0168781,-6.16236e-06,1.127e-09,-8.03422e-14,78980.3,-43.758], Tmin=(990.714,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(685.636,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(274.378,'J/(mol*K)'), comment="""Thermo group additivity estimation: group(Cs-(Cds-O2d)CsCsCs) + group(Cs-(Cds-Cds)CsHH) + group(Cs-CsHHH) + group(Cds-CdsCsCs) + group(Cds-OdCsH) + group(Cds-CdsHH) + ring(Methylene_cyclopropane) + radical(CCJ2_triplet) + radical(CCJCC=O)"""),
)

species(
    label = '[CH]C1[CH]OC(=C)C=1(24910)',
    structure = SMILES('[CH]C1[CH]OC(=C)C=1'),
    E0 = (398.371,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
    molecularWeight = (93.1033,'amu'),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[1.83264,0.0244272,7.70986e-05,-1.24228e-07,5.12932e-11,48012.3,20.1064], Tmin=(100,'K'), Tmax=(936.189,'K')), NASAPolynomial(coeffs=[17.7398,0.0155799,-3.44755e-06,5.81333e-10,-4.82628e-14,42443.2,-69.4238], Tmin=(936.189,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(398.371,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(278.535,'J/(mol*K)'), comment="""Thermo group additivity estimation: group(O2s-Cs(Cds-Cd)) + group(Cs-(Cds-Cds)OsHH) + group(Cs-(Cds-Cds)HHH) + group(Cds-CdsCsCs) + group(Cds-Cds(Cds-Cds)O2s) + group(Cds-Cds(Cds-Cds)H) + group(Cds-CdsHH) + ring(Cyclopentane) + radical(C=CCJ(O)C) + radical(AllylJ2_triplet)"""),
)

species(
    label = '[CH]C([C]=O)C=[C][CH2](24941)',
    structure = SMILES('[CH]C([C]=O)C=[C][CH2]'),
    E0 = (830.222,'kJ/mol'),
    modes = [
        HarmonicOscillator(frequencies=([3000,3100,440,815,1455,1000,1685,370,1380,1390,370,380,2900,435,1855,455,950,3010,987.5,1337.5,450,1655,180,794.886,3092.83,3092.83],'cm^-1')),
        HinderedRotor(inertia=(0.593924,'amu*angstrom^2'), symmetry=1, barrier=(13.6555,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.0993561,'amu*angstrom^2'), symmetry=1, barrier=(44.547,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.593937,'amu*angstrom^2'), symmetry=1, barrier=(13.6558,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.379483,'amu*angstrom^2'), symmetry=1, barrier=(44.547,'kJ/mol'), semiclassical=False),
    ],
    spinMultiplicity = 6,
    opticalIsomers = 1,
    molecularWeight = (93.1033,'amu'),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[1.21406,0.0652205,-7.99365e-05,5.53728e-08,-1.57867e-11,99949.5,28.281], Tmin=(100,'K'), Tmax=(847.785,'K')), NASAPolynomial(coeffs=[9.28234,0.0271532,-1.25838e-05,2.40946e-09,-1.68606e-13,98581.5,-9.31148], Tmin=(847.785,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(830.222,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(266.063,'J/(mol*K)'), comment="""Thermo group additivity estimation: group(Cs-(Cds-O2d)CsCsH) + group(Cs-CsHHH) + group(Cs-(Cds-Cds)HHH) + group(Cds-CdsCsH) + group(Cds-CdsCsH) + group(Cds-OdCsH) + radical(Allyl_P) + radical(CCJ2_triplet) + radical(CC(C)CJ=O) + radical(Cds_S)"""),
)

species(
    label = '[CH]C(=[C][C]=C)C[O](24942)',
    structure = SMILES('[CH]C(=[C][C]=C)C[O]'),
    E0 = (860.085,'kJ/mol'),
    modes = [
        HarmonicOscillator(frequencies=([1670,1700,300,440,2750,2850,1437.5,1250,1305,750,350,350,440,435,1725,2950,3100,1380,975,1025,1650,334.989,335.035,335.044,335.092,335.095,335.099],'cm^-1')),
        HinderedRotor(inertia=(0.663501,'amu*angstrom^2'), symmetry=1, barrier=(52.8525,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.663548,'amu*angstrom^2'), symmetry=1, barrier=(52.8528,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.663368,'amu*angstrom^2'), symmetry=1, barrier=(52.8525,'kJ/mol'), semiclassical=False),
    ],
    spinMultiplicity = 6,
    opticalIsomers = 1,
    molecularWeight = (93.1033,'amu'),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[0.869955,0.077066,-0.000113603,9.99245e-08,-3.43433e-11,103549,26.7112], Tmin=(100,'K'), Tmax=(886.845,'K')), NASAPolynomial(coeffs=[5.08892,0.0389794,-1.69502e-05,3.03711e-09,-2.00326e-13,103550,11.089], Tmin=(886.845,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(860.085,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(270.22,'J/(mol*K)'), comment="""Thermo group additivity estimation: group(O2s-CsH) + group(Cs-(Cds-Cds)OsHH) + group(Cs-(Cds-Cds)HHH) + group(Cds-CdsCsCs) + group(Cds-Cds(Cds-Cds)H) + group(Cds-Cds(Cds-Cds)H) + group(Cds-CdsHH) + radical(CCOJ) + radical(AllylJ2_triplet) + radical(C=CJC=C) + radical(C=CJC=C)"""),
)

species(
    label = '[CH]C([C]=[C][CH2])=CO(24943)',
    structure = SMILES('[CH]C([C]=[C][CH2])=CO'),
    E0 = (737.09,'kJ/mol'),
    modes = [
        HarmonicOscillator(frequencies=([3000,3100,440,815,1455,1000,1670,1700,300,440,3010,987.5,1337.5,450,1655,3615,1277.5,1000,350,440,435,1725,180,180,180,180],'cm^-1')),
        HinderedRotor(inertia=(2.09933,'amu*angstrom^2'), symmetry=1, barrier=(48.2678,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(2.09963,'amu*angstrom^2'), symmetry=1, barrier=(48.2747,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(2.10034,'amu*angstrom^2'), symmetry=1, barrier=(48.291,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(2.09928,'amu*angstrom^2'), symmetry=1, barrier=(48.2666,'kJ/mol'), semiclassical=False),
    ],
    spinMultiplicity = 6,
    opticalIsomers = 1,
    molecularWeight = (93.1033,'amu'),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[-0.516553,0.0803385,-8.31626e-05,4.24456e-08,-8.09632e-12,88830.3,29.2091], Tmin=(100,'K'), Tmax=(1497.05,'K')), NASAPolynomial(coeffs=[19.9602,0.0120882,-1.21319e-06,-8.87193e-11,1.54525e-14,84216.4,-72.7752], Tmin=(1497.05,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(737.09,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(266.063,'J/(mol*K)'), comment="""Thermo group additivity estimation: group(O2s-(Cds-Cd)H) + group(Cs-(Cds-Cds)HHH) + group(Cs-(Cds-Cds)HHH) + group(Cds-Cds(Cds-Cds)Cs) + group(Cds-CdsCsH) + group(Cds-Cds(Cds-Cds)H) + group(Cds-CdsOsH) + radical(Cds_S) + radical(AllylJ2_triplet) + radical(C=CJC=C) + radical(C=CC=CCJ)"""),
)

species(
    label = '[C-]#[O+](374)',
    structure = SMILES('[C-]#[O+]'),
    E0 = (299.89,'kJ/mol'),
    modes = [
        HarmonicOscillator(frequencies=([180],'cm^-1')),
    ],
    spinMultiplicity = 1,
    opticalIsomers = 1,
    molecularWeight = (28.0101,'amu'),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[3.33667,0.00896487,-2.66756e-05,3.61071e-08,-1.57199e-11,36069.2,-1.20266], Tmin=(100,'K'), Tmax=(865.594,'K')), NASAPolynomial(coeffs=[-0.394107,0.0117562,-6.47408e-06,1.26375e-09,-8.67562e-14,37256.3,19.3844], Tmin=(865.594,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(299.89,'kJ/mol'), Cp0=(29.1007,'J/(mol*K)'), CpInf=(37.4151,'J/(mol*K)'), comment="""Thermo group additivity estimation: group(O2s-CsCs) + group(CsJ2_singlet-CsH)"""),
)

species(
    label = '[CH]C=C[C]=C(17762)',
    structure = SMILES('[CH]C=C[C]=C'),
    E0 = (627.576,'kJ/mol'),
    modes = [
        HarmonicOscillator(frequencies=([2950,3100,1380,975,1025,1650,2995,3025,975,1000,1300,1375,400,500,1630,1680,1685,370,180,180,180,180],'cm^-1')),
        HinderedRotor(inertia=(2.14426,'amu*angstrom^2'), symmetry=1, barrier=(49.3007,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(2.12108,'amu*angstrom^2'), symmetry=1, barrier=(48.7679,'kJ/mol'), semiclassical=False),
    ],
    spinMultiplicity = 4,
    opticalIsomers = 1,
    molecularWeight = (65.0932,'amu'),
    collisionModel = TransportData(shapeIndex=2, epsilon=(3140.68,'J/mol'), sigma=(5.4037,'angstroms'), dipoleMoment=(0,'C*m'), polarizability=(0,'angstroms^3'), rotrelaxcollnum=0, comment="""Epsilon & sigma estimated with Tc=490.57 K, Pc=45.16 bar (from Joback method)"""),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[1.94063,0.0391755,-1.26068e-05,-1.10115e-08,7.27021e-12,75559.6,17.8981], Tmin=(100,'K'), Tmax=(954.651,'K')), NASAPolynomial(coeffs=[9.55667,0.0216437,-7.65371e-06,1.30771e-09,-8.77958e-14,73450.2,-21.9233], Tmin=(954.651,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(627.576,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(224.491,'J/(mol*K)'), comment="""Thermo group additivity estimation: group(Cs-(Cds-Cds)HHH) + group(Cds-CdsCsH) + group(Cds-Cds(Cds-Cds)H) + group(Cds-Cds(Cds-Cds)H) + group(Cds-CdsHH) + radical(AllylJ2_triplet) + radical(C=CJC=C)"""),
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
    label = '[CH]C(=[CH])C=O(21223)',
    structure = SMILES('[CH]C(=[CH])C=O'),
    E0 = (493.098,'kJ/mol'),
    modes = [
        HarmonicOscillator(frequencies=([2782.5,750,1395,475,1775,1000,350,440,435,1725,3120,650,792.5,1650,180,180],'cm^-1')),
        HinderedRotor(inertia=(2.15191,'amu*angstrom^2'), symmetry=1, barrier=(49.4767,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(2.15169,'amu*angstrom^2'), symmetry=1, barrier=(49.4716,'kJ/mol'), semiclassical=False),
    ],
    spinMultiplicity = 4,
    opticalIsomers = 1,
    molecularWeight = (67.066,'amu'),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[2.81344,0.0328631,-1.91319e-05,4.81753e-09,-4.76694e-13,59341.3,14.8318], Tmin=(100,'K'), Tmax=(2138.94,'K')), NASAPolynomial(coeffs=[10.3328,0.0188013,-9.27062e-06,1.74397e-09,-1.17457e-13,56124.6,-27.162], Tmin=(2138.94,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(493.098,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(174.604,'J/(mol*K)'), comment="""Thermo group additivity estimation: group(Cs-(Cds-Cds)HHH) + group(Cd-CdCs(CO)) + group(Cds-O2d(Cds-Cds)H) + group(Cds-CdsHH) + radical(AllylJ2_triplet) + radical(Cds_P)"""),
)

species(
    label = '[CH]C1=COC1[C]=C(24944)',
    structure = SMILES('[CH]C1=COC1[C]=C'),
    E0 = (632.219,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
    molecularWeight = (93.1033,'amu'),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[1.0574,0.0431634,3.30898e-05,-8.60652e-08,3.9302e-11,76164.1,23.2338], Tmin=(100,'K'), Tmax=(940.322,'K')), NASAPolynomial(coeffs=[21.4676,0.010751,-2.00059e-06,3.48376e-10,-3.26443e-14,69920.2,-86.7684], Tmin=(940.322,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(632.219,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(274.378,'J/(mol*K)'), comment="""Thermo group additivity estimation: group(O2s-Cs(Cds-Cd)) + group(Cs-(Cds-Cds)(Cds-Cds)OsH) + group(Cs-(Cds-Cds)HHH) + group(Cds-CdsCsCs) + group(Cds-CdsCsH) + group(Cds-CdsOsH) + group(Cds-CdsHH) + ring(Cyclobutene) + radical(AllylJ2_triplet) + radical(Cds_S)"""),
)

species(
    label = '[CH]C(=C=O)C=C=C(24945)',
    structure = SMILES('[CH]C(=C=O)C=C=C'),
    E0 = (511.333,'kJ/mol'),
    modes = [
        HarmonicOscillator(frequencies=([2950,3100,1380,975,1025,1650,540,610,2055,350,440,435,1725,2120,512.5,787.5,3010,987.5,1337.5,450,1655,180,180,180,180],'cm^-1')),
        HinderedRotor(inertia=(2.18362,'amu*angstrom^2'), symmetry=1, barrier=(50.2058,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(2.18454,'amu*angstrom^2'), symmetry=1, barrier=(50.227,'kJ/mol'), semiclassical=False),
    ],
    spinMultiplicity = 3,
    opticalIsomers = 1,
    molecularWeight = (92.0954,'amu'),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[2.03517,0.0522042,-1.87105e-05,-6.91697e-08,7.90623e-11,61560.6,20.5556], Tmin=(100,'K'), Tmax=(466.34,'K')), NASAPolynomial(coeffs=[5.95313,0.0350397,-1.63851e-05,3.10843e-09,-2.14942e-13,61016.4,2.72567], Tmin=(466.34,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(511.333,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(249.434,'J/(mol*K)'), comment="""Thermo group additivity estimation: group(Cs-(Cds-Cdd-O2d)HHH) + group(Cds-(Cdd-O2d)CsCs) + group(Cd-Cd(CCO)H) + group(Cds-CdsHH) + group(Cdd-CdsCds) + radical(AllylJ2_triplet)"""),
)

species(
    label = '[CH]C(=[C]O)C=C=C(24946)',
    structure = SMILES('[CH]C(=[C]O)C=C=C'),
    E0 = (598.545,'kJ/mol'),
    modes = [
        HarmonicOscillator(frequencies=([3615,1277.5,1000,540,610,2055,350,440,435,1725,2950,3100,1380,975,1025,1650,1685,370,3010,987.5,1337.5,450,1655,180,180,180,180],'cm^-1')),
        HinderedRotor(inertia=(2.00077,'amu*angstrom^2'), symmetry=1, barrier=(46.0017,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(1.99184,'amu*angstrom^2'), symmetry=1, barrier=(45.7964,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(1.99872,'amu*angstrom^2'), symmetry=1, barrier=(45.9545,'kJ/mol'), semiclassical=False),
    ],
    spinMultiplicity = 4,
    opticalIsomers = 1,
    molecularWeight = (93.1033,'amu'),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[-0.189093,0.0765065,-7.76001e-05,3.87268e-08,-7.33033e-12,72152.5,28.5926], Tmin=(100,'K'), Tmax=(1449.28,'K')), NASAPolynomial(coeffs=[19.9824,0.0124785,-2.68409e-06,2.8787e-10,-1.34778e-14,67183.2,-73.1809], Tmin=(1449.28,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(598.545,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(270.22,'J/(mol*K)'), comment="""Thermo group additivity estimation: group(O2s-(Cds-Cd)H) + group(Cs-(Cds-Cds)HHH) + group(Cds-Cds(Cds-Cds)Cs) + group(Cds-Cds(Cds-Cds)H) + group(Cds-CdsOsH) + group(Cds-CdsHH) + group(Cdd-CdsCds) + radical(AllylJ2_triplet) + radical(C=CJO)"""),
)

species(
    label = '[CH]C(C#C[CH2])=CO(24947)',
    structure = SMILES('[CH]C(C#C[CH2])=CO'),
    E0 = (506.369,'kJ/mol'),
    modes = [
        HarmonicOscillator(frequencies=([3000,3100,440,815,1455,1000,3010,987.5,1337.5,450,1655,3615,1277.5,1000,2100,2250,500,550,350,440,435,1725,180,180,180,180],'cm^-1')),
        HinderedRotor(inertia=(2.11603,'amu*angstrom^2'), symmetry=1, barrier=(48.6518,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(2.11869,'amu*angstrom^2'), symmetry=1, barrier=(48.7128,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(2.1132,'amu*angstrom^2'), symmetry=1, barrier=(48.5865,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(2.11981,'amu*angstrom^2'), symmetry=1, barrier=(48.7385,'kJ/mol'), semiclassical=False),
    ],
    spinMultiplicity = 4,
    opticalIsomers = 1,
    molecularWeight = (93.1033,'amu'),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[0.765853,0.0584163,-2.1872e-05,-2.46111e-08,1.7135e-11,61030.3,25.5165], Tmin=(100,'K'), Tmax=(924.64,'K')), NASAPolynomial(coeffs=[18.5533,0.0142791,-3.49861e-06,5.19375e-10,-3.60631e-14,56338.3,-66.4885], Tmin=(924.64,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(506.369,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(266.063,'J/(mol*K)'), comment="""Thermo group additivity estimation: group(O2s-(Cds-Cd)H) + group(Cs-(Cds-Cds)HHH) + group(Cs-CtHHH) + group(Cds-CdsCtCs) + group(Cds-CdsOsH) + group(Ct-CtCs) + group(Ct-Ct(Cds-Cds)) + radical(Propargyl) + radical(AllylJ2_triplet)"""),
)

species(
    label = '[CH]C(C=C=[CH])=CO(24948)',
    structure = SMILES('[CH]C(C=C=[CH])=CO'),
    E0 = (513.278,'kJ/mol'),
    modes = [
        HarmonicOscillator(frequencies=([3120,650,792.5,1650,540,610,2055,2995,3025,975,1000,1300,1375,400,500,1630,1680,3615,1277.5,1000,350,440,435,1725,180,180,180],'cm^-1')),
        HinderedRotor(inertia=(1.96074,'amu*angstrom^2'), symmetry=1, barrier=(45.0813,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(1.96208,'amu*angstrom^2'), symmetry=1, barrier=(45.1121,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(1.96077,'amu*angstrom^2'), symmetry=1, barrier=(45.0819,'kJ/mol'), semiclassical=False),
    ],
    spinMultiplicity = 4,
    opticalIsomers = 1,
    molecularWeight = (93.1033,'amu'),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[-1.12477,0.0843282,-8.82067e-05,4.37212e-08,-7.92205e-12,61942.1,29.4583], Tmin=(100,'K'), Tmax=(1625.88,'K')), NASAPolynomial(coeffs=[22.2226,0.00699829,1.48662e-06,-5.8048e-10,4.6859e-14,56979.1,-86.4422], Tmin=(1625.88,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(513.278,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(270.22,'J/(mol*K)'), comment="""Thermo group additivity estimation: group(O2s-(Cds-Cd)H) + group(Cs-(Cds-Cds)HHH) + group(Cds-Cds(Cds-Cds)Cs) + group(Cds-Cds(Cds-Cds)H) + group(Cds-CdsOsH) + group(Cds-CdsHH) + group(Cdd-CdsCds) + radical(AllylJ2_triplet) + radical(C=C=CJ)"""),
)

species(
    label = '[CH][C]=C[O](21209)',
    structure = SMILES('[CH][C]=C[O]'),
    E0 = (547.32,'kJ/mol'),
    modes = [
        HarmonicOscillator(frequencies=([1685,370,3010,987.5,1337.5,450,1655,180,180,180,180],'cm^-1')),
        HinderedRotor(inertia=(2.12843,'amu*angstrom^2'), symmetry=1, barrier=(48.9368,'kJ/mol'), semiclassical=False),
    ],
    spinMultiplicity = 5,
    opticalIsomers = 1,
    molecularWeight = (54.0474,'amu'),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[2.67099,0.0213545,9.09852e-06,-3.1272e-08,1.4879e-11,65882.5,14.7882], Tmin=(100,'K'), Tmax=(925.361,'K')), NASAPolynomial(coeffs=[10.425,0.00802333,-2.01432e-06,3.08705e-10,-2.20542e-14,63583.2,-26.689], Tmin=(925.361,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(547.32,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(128.874,'J/(mol*K)'), comment="""Thermo group additivity estimation: group(O2s-(Cds-Cd)H) + group(Cs-(Cds-Cds)HHH) + group(Cds-CdsCsH) + group(Cds-CdsOsH) + radical(C=COJ) + radical(Cds_S) + radical(AllylJ2_triplet)"""),
)

species(
    label = '[CH]C1([CH]O1)C=C=C(24949)',
    structure = SMILES('[CH]C1([CH]O1)C=C=C'),
    E0 = (723.931,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
    molecularWeight = (93.1033,'amu'),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[-0.36234,0.0728903,-7.80113e-05,4.08047e-08,-7.75216e-12,87245.8,28.6338], Tmin=(100,'K'), Tmax=(1579.63,'K')), NASAPolynomial(coeffs=[17.4698,0.00822382,1.92283e-06,-7.49922e-10,6.1995e-14,84046.4,-57.8437], Tmin=(1579.63,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(723.931,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(274.378,'J/(mol*K)'), comment="""Thermo group additivity estimation: group(O2s-CsCs) + group(Cs-(Cds-Cds)CsCsOs) + group(Cs-CsOsHH) + group(Cs-CsHHH) + group(Cds-CdsCsH) + group(Cds-CdsHH) + group(Cdd-CdsCds) + ring(Ethylene_oxide) + radical(CCJ2_triplet) + radical(CCsJO)"""),
)

species(
    label = '[CH]C1C=[C]COC=1(24822)',
    structure = SMILES('[CH]C1C=[C]COC=1'),
    E0 = (519.66,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
    molecularWeight = (93.1033,'amu'),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[1.34816,0.037668,4.27034e-05,-9.17473e-08,4.05137e-11,62615.1,18.4099], Tmin=(100,'K'), Tmax=(936.683,'K')), NASAPolynomial(coeffs=[19.1915,0.0138488,-3.03159e-06,5.03072e-10,-4.15211e-14,56974.6,-78.7724], Tmin=(936.683,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(519.66,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(278.535,'J/(mol*K)'), comment="""Thermo group additivity estimation: group(O2s-Cs(Cds-Cd)) + group(Cs-(Cds-Cds)OsHH) + group(Cs-(Cds-Cds)HHH) + group(Cds-Cds(Cds-Cds)Cs) + group(Cds-CdsCsH) + group(Cds-Cds(Cds-Cds)H) + group(Cds-CdsOsH) + ring(1,3-Cyclohexadiene) + radical(AllylJ2_triplet) + radical(Cds_S)"""),
)

species(
    label = '[CH]C(C#C[CH2])[CH][O](24950)',
    structure = SMILES('[CH]C(C#C[CH2])[CH][O]'),
    E0 = (925.293,'kJ/mol'),
    modes = [
        HarmonicOscillator(frequencies=([3025,407.5,1350,352.5,1380,1390,370,380,2900,435,2100,2250,500,550,3000,3100,440,815,1455,1000,200,800,1000,1200,1400,1600],'cm^-1')),
        HinderedRotor(inertia=(0.156089,'amu*angstrom^2'), symmetry=1, barrier=(3.5888,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.156089,'amu*angstrom^2'), symmetry=1, barrier=(3.5888,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.156089,'amu*angstrom^2'), symmetry=1, barrier=(3.5888,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.156089,'amu*angstrom^2'), symmetry=1, barrier=(3.5888,'kJ/mol'), semiclassical=False),
    ],
    spinMultiplicity = 6,
    opticalIsomers = 1,
    molecularWeight = (93.1033,'amu'),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[1.02726,0.0714976,-0.000106364,9.01692e-08,-3.01147e-11,111388,29.151], Tmin=(100,'K'), Tmax=(866.955,'K')), NASAPolynomial(coeffs=[7.45315,0.029743,-1.31733e-05,2.40078e-09,-1.60478e-13,110729,1.69114], Tmin=(866.955,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(925.293,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(266.063,'J/(mol*K)'), comment="""Thermo group additivity estimation: group(O2s-CsH) + group(Cs-CtCsCsH) + group(Cs-CsOsHH) + group(Cs-CsHHH) + group(Cs-CtHHH) + group(Ct-CtCs) + group(Ct-CtCs) + radical(CCsJOH) + radical(Propargyl) + radical(CCOJ) + radical(CCJ2_triplet)"""),
)

species(
    label = '[CH]C([C]=C[CH2])=C[O](24951)',
    structure = SMILES('[CH]C([C]=C[CH2])=C[O]'),
    E0 = (640.711,'kJ/mol'),
    modes = [
        HarmonicOscillator(frequencies=([3000,3100,440,815,1455,1000,1685,370,350,440,435,1725,2995,3025,975,1000,1300,1375,400,500,1630,1680,180,180,180,180,180],'cm^-1')),
        HinderedRotor(inertia=(2.14828,'amu*angstrom^2'), symmetry=1, barrier=(49.3931,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(2.14822,'amu*angstrom^2'), symmetry=1, barrier=(49.3918,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(2.1482,'amu*angstrom^2'), symmetry=1, barrier=(49.3914,'kJ/mol'), semiclassical=False),
    ],
    spinMultiplicity = 6,
    opticalIsomers = 1,
    molecularWeight = (93.1033,'amu'),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[0.78553,0.0586322,-1.93178e-05,-2.73233e-08,1.83398e-11,77186.7,25.474], Tmin=(100,'K'), Tmax=(910.036,'K')), NASAPolynomial(coeffs=[17.5645,0.0171698,-4.19617e-06,5.86769e-10,-3.80161e-14,72795.8,-61.2388], Tmin=(910.036,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(640.711,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(270.22,'J/(mol*K)'), comment="""Thermo group additivity estimation: group(O2s-(Cds-Cd)H) + group(Cs-(Cds-Cds)HHH) + group(Cs-(Cds-Cds)HHH) + group(Cds-Cds(Cds-Cds)Cs) + group(Cds-CdsCsH) + group(Cds-Cds(Cds-Cds)H) + group(Cds-CdsOsH) + radical(C=COJ) + radical(AllylJ2_triplet) + radical(C=CJC=C) + radical(C=CC=CCJ)"""),
)

species(
    label = '[CH]C=CC([CH])=C[O](24952)',
    structure = SMILES('[CH]C=CC([CH])=C[O]'),
    E0 = (694.344,'kJ/mol'),
    modes = [
        HarmonicOscillator(frequencies=([350,440,435,1725,2995,3010,3025,975,987.5,1000,1300,1337.5,1375,400,450,500,1630,1655,1680,200,800,933.333,1066.67,1200,1333.33,1466.67,1600],'cm^-1')),
        HinderedRotor(inertia=(0.156089,'amu*angstrom^2'), symmetry=1, barrier=(3.5888,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.156089,'amu*angstrom^2'), symmetry=1, barrier=(3.5888,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.156089,'amu*angstrom^2'), symmetry=1, barrier=(3.5888,'kJ/mol'), semiclassical=False),
    ],
    spinMultiplicity = 6,
    opticalIsomers = 1,
    molecularWeight = (93.1033,'amu'),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[0.576379,0.0610234,-1.36041e-05,-3.39224e-08,1.98315e-11,83646.6,26.1773], Tmin=(100,'K'), Tmax=(941.254,'K')), NASAPolynomial(coeffs=[18.2189,0.0213974,-6.78704e-06,1.14751e-09,-8.02837e-14,78759.5,-66.1873], Tmin=(941.254,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(694.344,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(270.22,'J/(mol*K)'), comment="""Thermo group additivity estimation: group(O2s-(Cds-Cd)H) + group(Cs-(Cds-Cds)HHH) + group(Cs-(Cds-Cds)HHH) + group(Cds-Cds(Cds-Cds)Cs) + group(Cds-CdsCsH) + group(Cds-Cds(Cds-Cds)H) + group(Cds-CdsOsH) + radical(AllylJ2_triplet) + radical(AllylJ2_triplet) + radical(C=COJ)"""),
)

species(
    label = '[CH]C(=[C][O])C[C]=C(24953)',
    structure = SMILES('[CH]C(=[C][O])C[C]=C'),
    E0 = (830.635,'kJ/mol'),
    modes = [
        HarmonicOscillator(frequencies=([1670,1700,300,440,2750,2850,1437.5,1250,1305,750,350,350,440,435,1725,2950,3100,1380,975,1025,1650,410.096,410.096,410.099,410.103,410.104,410.106],'cm^-1')),
        HinderedRotor(inertia=(0.429195,'amu*angstrom^2'), symmetry=1, barrier=(51.2233,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.429208,'amu*angstrom^2'), symmetry=1, barrier=(51.2233,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.429203,'amu*angstrom^2'), symmetry=1, barrier=(51.2233,'kJ/mol'), semiclassical=False),
    ],
    spinMultiplicity = 6,
    opticalIsomers = 1,
    molecularWeight = (93.1033,'amu'),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[1.08946,0.0614052,-5.28599e-05,2.4461e-08,-4.66053e-12,100009,30.0083], Tmin=(100,'K'), Tmax=(1240.7,'K')), NASAPolynomial(coeffs=[11.8078,0.0268491,-1.10815e-05,2.01206e-09,-1.37049e-13,97349.6,-24.0133], Tmin=(1240.7,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(830.635,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(270.22,'J/(mol*K)'), comment="""Thermo group additivity estimation: group(O2s-(Cds-Cd)H) + group(Cs-(Cds-Cds)(Cds-Cds)HH) + group(Cs-(Cds-Cds)HHH) + group(Cds-CdsCsCs) + group(Cds-CdsCsH) + group(Cds-CdsOsH) + group(Cds-CdsHH) + radical(AllylJ2_triplet) + radical(Cds_S) + radical(C=CJO) + radical(C=COJ)"""),
)

species(
    label = '[CH]C(=[C][O])C=C[CH2](24954)',
    structure = SMILES('[CH]C(=[C][O])C=C[CH2]'),
    E0 = (681.459,'kJ/mol'),
    modes = [
        HarmonicOscillator(frequencies=([3000,3100,440,815,1455,1000,1685,370,2995,3025,975,1000,1300,1375,400,500,1630,1680,350,440,435,1725,313.672,313.673,313.673,313.674,313.674],'cm^-1')),
        HinderedRotor(inertia=(0.722172,'amu*angstrom^2'), symmetry=1, barrier=(50.4219,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.722166,'amu*angstrom^2'), symmetry=1, barrier=(50.4219,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.722171,'amu*angstrom^2'), symmetry=1, barrier=(50.422,'kJ/mol'), semiclassical=False),
    ],
    spinMultiplicity = 6,
    opticalIsomers = 1,
    molecularWeight = (93.1033,'amu'),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[1.0578,0.054055,-1.79675e-05,-2.00723e-08,1.31933e-11,82076.3,28.0271], Tmin=(100,'K'), Tmax=(946.267,'K')), NASAPolynomial(coeffs=[14.9708,0.0211943,-7.01493e-06,1.19368e-09,-8.21335e-14,78281.3,-44.4659], Tmin=(946.267,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(681.459,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(270.22,'J/(mol*K)'), comment="""Thermo group additivity estimation: group(O2s-(Cds-Cd)H) + group(Cs-(Cds-Cds)HHH) + group(Cs-(Cds-Cds)HHH) + group(Cds-Cds(Cds-Cds)Cs) + group(Cds-CdsCsH) + group(Cds-Cds(Cds-Cds)H) + group(Cds-CdsOsH) + radical(C=COJ) + radical(C=CC=CCJ) + radical(AllylJ2_triplet) + radical(C=CJO)"""),
)

species(
    label = '[CH][C]=CC([CH])C=O(21904)',
    structure = SMILES('[CH][C]=CC([CH])C=O'),
    E0 = (890.719,'kJ/mol'),
    modes = [
        HarmonicOscillator(frequencies=([1685,370,2782.5,750,1395,475,1775,1000,1380,1390,370,380,2900,435,3010,987.5,1337.5,450,1655,200,800,960,1120,1280,1440,1600],'cm^-1')),
        HinderedRotor(inertia=(0.156089,'amu*angstrom^2'), symmetry=1, barrier=(3.5888,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.156089,'amu*angstrom^2'), symmetry=1, barrier=(3.5888,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.156089,'amu*angstrom^2'), symmetry=1, barrier=(3.5888,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.156089,'amu*angstrom^2'), symmetry=1, barrier=(3.5888,'kJ/mol'), semiclassical=False),
    ],
    spinMultiplicity = 6,
    opticalIsomers = 1,
    molecularWeight = (93.1033,'amu'),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[1.04597,0.070439,-9.21982e-05,7.6427e-08,-2.65117e-11,107230,28.5458], Tmin=(100,'K'), Tmax=(795.179,'K')), NASAPolynomial(coeffs=[6.28358,0.0374637,-1.74907e-05,3.31036e-09,-2.2838e-13,106606,5.79547], Tmin=(795.179,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(890.719,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(266.063,'J/(mol*K)'), comment="""Thermo group additivity estimation: group(Cs-(Cds-O2d)CsCsH) + group(Cs-CsHHH) + group(Cs-(Cds-Cds)HHH) + group(Cds-CdsCsH) + group(Cds-CdsCsH) + group(Cds-OdCsH) + radical(AllylJ2_triplet) + radical(CCJ2_triplet) + radical(Cds_S)"""),
)

species(
    label = '[CH][C]([C]=O)C=[C]C(24955)',
    structure = SMILES('[CH][C]([C]=O)C=[C]C'),
    E0 = (798.424,'kJ/mol'),
    modes = [
        HarmonicOscillator(frequencies=([1685,370,2750,2800,2850,1350,1500,750,1050,1375,1000,3010,987.5,1337.5,450,1655,1855,455,950,360,370,350,397.745,397.753,397.852,397.949],'cm^-1')),
        HinderedRotor(inertia=(0.00106476,'amu*angstrom^2'), symmetry=1, barrier=(0.119627,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.0010661,'amu*angstrom^2'), symmetry=1, barrier=(0.119627,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.00106495,'amu*angstrom^2'), symmetry=1, barrier=(0.119627,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.00106545,'amu*angstrom^2'), symmetry=1, barrier=(0.119627,'kJ/mol'), semiclassical=False),
    ],
    spinMultiplicity = 6,
    opticalIsomers = 1,
    molecularWeight = (93.1033,'amu'),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[1.69051,0.0512712,-4.14559e-05,1.79255e-08,-3.26089e-12,96110.8,25.7363], Tmin=(100,'K'), Tmax=(1265.49,'K')), NASAPolynomial(coeffs=[9.66225,0.026074,-1.15896e-05,2.19185e-09,-1.52695e-13,94093.2,-14.5998], Tmin=(1265.49,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(798.424,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(266.063,'J/(mol*K)'), comment="""Thermo group additivity estimation: group(Cs-(Cds-O2d)CsCsH) + group(Cs-CsHHH) + group(Cs-(Cds-Cds)HHH) + group(Cds-CdsCsH) + group(Cds-CdsCsH) + group(Cds-OdCsH) + radical(Cds_S) + radical(CCJ2_triplet) + radical(CC(C)CJ=O) + radical(C=CCJ(C)C=O)"""),
)

species(
    label = '[CH][C]([CH]C#C)C[O](24779)',
    structure = SMILES('[CH][C]([CH]C#C)C[O]'),
    E0 = (905.857,'kJ/mol'),
    modes = [
        HarmonicOscillator(frequencies=([3025,407.5,1350,352.5,2750,2850,1437.5,1250,1305,750,350,750,770,3400,2100,2175,525,360,370,350,200,800,1000,1200,1400,1600],'cm^-1')),
        HinderedRotor(inertia=(0.156089,'amu*angstrom^2'), symmetry=1, barrier=(3.5888,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.156089,'amu*angstrom^2'), symmetry=1, barrier=(3.5888,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.156089,'amu*angstrom^2'), symmetry=1, barrier=(3.5888,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.156089,'amu*angstrom^2'), symmetry=1, barrier=(3.5888,'kJ/mol'), semiclassical=False),
    ],
    spinMultiplicity = 6,
    opticalIsomers = 1,
    molecularWeight = (93.1033,'amu'),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[1.1978,0.0651071,-8.52633e-05,6.60118e-08,-2.05617e-11,109047,26.7905], Tmin=(100,'K'), Tmax=(888.202,'K')), NASAPolynomial(coeffs=[7.90985,0.0281201,-1.13841e-05,1.99133e-09,-1.30351e-13,108121,-3.2946], Tmin=(888.202,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(905.857,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(266.063,'J/(mol*K)'), comment="""Thermo group additivity estimation: group(O2s-CsH) + group(Cs-CsCsCsH) + group(Cs-CtCsHH) + group(Cs-CsOsHH) + group(Cs-CsHHH) + group(Ct-CtCs) + group(Ct-CtH) + radical(CCJ(C)CO) + radical(CCJ2_triplet) + radical(CCOJ) + radical(Sec_Propargyl)"""),
)

species(
    label = '[CH]C(=CC#C)C[O](24956)',
    structure = SMILES('[CH]C(=CC#C)C[O]'),
    E0 = (638.88,'kJ/mol'),
    modes = [
        HarmonicOscillator(frequencies=([2750,2850,1437.5,1250,1305,750,350,750,770,3400,2100,350,440,435,1725,2175,525,3010,987.5,1337.5,450,1655,369.123,369.124,369.124,369.124,369.125],'cm^-1')),
        HinderedRotor(inertia=(0.530029,'amu*angstrom^2'), symmetry=1, barrier=(51.2473,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.53003,'amu*angstrom^2'), symmetry=1, barrier=(51.2473,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.53003,'amu*angstrom^2'), symmetry=1, barrier=(51.2473,'kJ/mol'), semiclassical=False),
    ],
    spinMultiplicity = 4,
    opticalIsomers = 1,
    molecularWeight = (93.1033,'amu'),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[1.50502,0.0594049,-5.51419e-05,3.26869e-08,-8.85115e-12,76925.5,24.1494], Tmin=(100,'K'), Tmax=(849.2,'K')), NASAPolynomial(coeffs=[5.98259,0.0383148,-1.78902e-05,3.44328e-09,-2.42269e-13,76165,3.27951], Tmin=(849.2,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(638.88,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(270.22,'J/(mol*K)'), comment="""Thermo group additivity estimation: group(O2s-CsH) + group(Cs-(Cds-Cds)OsHH) + group(Cs-(Cds-Cds)HHH) + group(Cds-CdsCsCs) + group(Cds-CdsCtH) + group(Ct-Ct(Cds-Cds)) + group(Ct-CtH) + radical(AllylJ2_triplet) + radical(CCOJ)"""),
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
    label = '[CH]C(=[CH])C=C=C(21246)',
    structure = SMILES('[CH]C(=[CH])C=C=C'),
    E0 = (814.69,'kJ/mol'),
    modes = [
        HarmonicOscillator(frequencies=([3120,650,792.5,1650,2950,3100,1380,975,1025,1650,540,610,2055,350,440,435,1725,3010,987.5,1337.5,450,1655,180,180,180],'cm^-1')),
        HinderedRotor(inertia=(2.04822,'amu*angstrom^2'), symmetry=1, barrier=(47.0927,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(2.0482,'amu*angstrom^2'), symmetry=1, barrier=(47.0921,'kJ/mol'), semiclassical=False),
    ],
    spinMultiplicity = 4,
    opticalIsomers = 1,
    molecularWeight = (77.1039,'amu'),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[1.13547,0.0558545,-3.99179e-05,9.1035e-09,1.39521e-12,98094,20.5557], Tmin=(100,'K'), Tmax=(1002.55,'K')), NASAPolynomial(coeffs=[13.4638,0.0204617,-7.60395e-06,1.34079e-09,-9.16303e-14,94928.7,-42.4104], Tmin=(1002.55,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(814.69,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(249.434,'J/(mol*K)'), comment="""Thermo group additivity estimation: group(Cs-(Cds-Cds)HHH) + group(Cds-Cds(Cds-Cds)Cs) + group(Cds-Cds(Cds-Cds)H) + group(Cds-CdsHH) + group(Cds-CdsHH) + group(Cdd-CdsCds) + radical(Cds_P) + radical(AllylJ2_triplet)"""),
)

species(
    label = 'C=[C]C=C1[CH]C1[O](24957)',
    structure = SMILES('C=[C]C=C1[CH]C1[O]'),
    E0 = (612.183,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
    molecularWeight = (93.1033,'amu'),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[0.922226,0.0530055,-1.13055e-05,-3.51083e-08,2.05865e-11,73752.8,21.0365], Tmin=(100,'K'), Tmax=(938.386,'K')), NASAPolynomial(coeffs=[19.5734,0.0104126,-2.22224e-06,3.55182e-10,-2.86535e-14,68627.3,-76.4178], Tmin=(938.386,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(612.183,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(278.535,'J/(mol*K)'), comment="""Thermo group additivity estimation: group(O2s-CsH) + group(Cs-(Cds-Cds)CsOsH) + group(Cs-(Cds-Cds)CsHH) + group(Cds-CdsCsCs) + group(Cds-Cds(Cds-Cds)H) + group(Cds-Cds(Cds-Cds)H) + group(Cds-CdsHH) + ring(Methylene_cyclopropane) + radical(C=CJC=C) + radical(CC(C)OJ) + radical(C=CCJCO)"""),
)

species(
    label = '[CH]=C1C=[C]CC1[O](24833)',
    structure = SMILES('[CH]=C1C=[C]CC1[O]'),
    E0 = (659.114,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
    molecularWeight = (93.1033,'amu'),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[1.68864,0.0358454,1.95116e-05,-5.50346e-08,2.44633e-11,79369.8,22.7905], Tmin=(100,'K'), Tmax=(968.517,'K')), NASAPolynomial(coeffs=[15.5596,0.0151743,-5.18391e-06,9.99888e-10,-7.67977e-14,74965.6,-52.5512], Tmin=(968.517,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(659.114,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(282.692,'J/(mol*K)'), comment="""Thermo group additivity estimation: group(O2s-CsH) + group(Cs-(Cds-Cds)CsOsH) + group(Cs-(Cds-Cds)CsHH) + group(Cds-Cds(Cds-Cds)Cs) + group(Cds-CdsCsH) + group(Cds-Cds(Cds-Cds)H) + group(Cds-CdsHH) + ring(3-Methylenecyclopentene) + radical(cyclopentene-vinyl) + radical(Cds_P) + radical(CC(C)OJ)"""),
)

species(
    label = '[CH][C]=C(18825)',
    structure = SMILES('[CH][C]=C'),
    E0 = (614.65,'kJ/mol'),
    modes = [
        HarmonicOscillator(frequencies=([2950,3100,1380,975,1025,1650,1685,370,228.264,228.889,229.07],'cm^-1')),
        HinderedRotor(inertia=(1.35219,'amu*angstrom^2'), symmetry=1, barrier=(50.6528,'kJ/mol'), semiclassical=False),
    ],
    spinMultiplicity = 4,
    opticalIsomers = 1,
    molecularWeight = (39.0559,'amu'),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[3.27541,0.0127954,9.49515e-06,-1.56026e-08,5.42938e-12,73954,11.3502], Tmin=(100,'K'), Tmax=(1063.31,'K')), NASAPolynomial(coeffs=[4.18965,0.0168435,-6.77763e-06,1.22218e-09,-8.33556e-14,73336.3,4.89309], Tmin=(1063.31,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(614.65,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(128.874,'J/(mol*K)'), comment="""Thermo group additivity estimation: group(Cs-(Cds-Cds)HHH) + group(Cds-CdsCsH) + group(Cds-CdsHH) + radical(Cds_S) + radical(AllylJ2_triplet)"""),
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
    label = '[CH2]C#CC([CH2])=C[O](24958)',
    structure = SMILES('[CH2]C#CC([CH2])=C[O]'),
    E0 = (428.647,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
    molecularWeight = (93.1033,'amu'),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[1.17091,0.0491069,-1.06115e-05,-3.07211e-08,1.79181e-11,51668.1,24.6441], Tmin=(100,'K'), Tmax=(943.907,'K')), NASAPolynomial(coeffs=[17.659,0.0118058,-3.09353e-06,5.25496e-10,-3.99664e-14,47104.5,-61.6352], Tmin=(943.907,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(428.647,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(270.22,'J/(mol*K)'), comment="""Thermo group additivity estimation: group(O2s-(Cds-Cd)H) + group(Cs-(Cds-Cds)HHH) + group(Cs-CtHHH) + group(Cds-CdsCtCs) + group(Cds-CdsOsH) + group(Ct-CtCs) + group(Ct-Ct(Cds-Cds)) + radical(Allyl_P) + radical(C=COJ) + radical(Propargyl)"""),
)

species(
    label = '[CH2]C([C]=O)=C[C]=C(24959)',
    structure = SMILES('[CH2]C([C]=O)=C[C]=C'),
    E0 = (446.041,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
    molecularWeight = (93.1033,'amu'),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[1.01321,0.0660113,-7.08233e-05,3.94135e-08,-8.79909e-12,53753.6,22.8095], Tmin=(100,'K'), Tmax=(1082.34,'K')), NASAPolynomial(coeffs=[12.9049,0.0220633,-9.91652e-06,1.89808e-09,-1.33766e-13,51179.5,-35.5021], Tmin=(1082.34,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(446.041,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(270.22,'J/(mol*K)'), comment="""Thermo group additivity estimation: group(Cs-(Cds-Cds)HHH) + group(Cd-CdCs(CO)) + group(Cds-Cds(Cds-Cds)H) + group(Cds-Cds(Cds-Cds)H) + group(Cds-O2d(Cds-Cds)H) + group(Cds-CdsHH) + radical(C=CJC=C) + radical(C=C(C)CJ=O) + radical(C=C(C=O)CJ)"""),
)

species(
    label = '[CH]C(C#CC)=C[O](24960)',
    structure = SMILES('[CH]C(C#CC)=C[O]'),
    E0 = (491.425,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
    molecularWeight = (93.1033,'amu'),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[1.16112,0.051508,-1.15786e-05,-2.79872e-08,1.67319e-11,59217.1,24.4518], Tmin=(100,'K'), Tmax=(922.978,'K')), NASAPolynomial(coeffs=[14.9989,0.0197343,-5.76523e-06,9.1274e-10,-6.12733e-14,55461.6,-47.7046], Tmin=(922.978,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(491.425,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(270.22,'J/(mol*K)'), comment="""Thermo group additivity estimation: group(O2s-(Cds-Cd)H) + group(Cs-(Cds-Cds)HHH) + group(Cs-CtHHH) + group(Cds-CdsCtCs) + group(Cds-CdsOsH) + group(Ct-CtCs) + group(Ct-Ct(Cds-Cds)) + radical(AllylJ2_triplet) + radical(C=COJ)"""),
)

species(
    label = '[CH]C(=C=O)C=[C]C(24961)',
    structure = SMILES('[CH]C(=C=O)C=[C]C'),
    E0 = (572.571,'kJ/mol'),
    modes = [
        HarmonicOscillator(frequencies=([1685,370,2750,2800,2850,1350,1500,750,1050,1375,1000,2120,512.5,787.5,350,440,435,1725,3010,987.5,1337.5,450,1655,344.058,346.368,347.847,350.594],'cm^-1')),
        HinderedRotor(inertia=(0.600746,'amu*angstrom^2'), symmetry=1, barrier=(53.2057,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.633162,'amu*angstrom^2'), symmetry=1, barrier=(53.2046,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.606508,'amu*angstrom^2'), symmetry=1, barrier=(53.1849,'kJ/mol'), semiclassical=False),
    ],
    spinMultiplicity = 4,
    opticalIsomers = 1,
    molecularWeight = (93.1033,'amu'),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[1.20584,0.0693282,-9.93984e-05,9.10896e-08,-3.35287e-11,68957.4,24.781], Tmin=(100,'K'), Tmax=(832.24,'K')), NASAPolynomial(coeffs=[3.70572,0.0411348,-1.94247e-05,3.66879e-09,-2.51505e-13,69101.5,16.5457], Tmin=(832.24,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(572.571,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(270.22,'J/(mol*K)'), comment="""Thermo group additivity estimation: group(Cs-(Cds-Cdd-O2d)HHH) + group(Cs-(Cds-Cds)HHH) + group(Cds-(Cdd-O2d)CsCs) + group(Cds-CdsCsH) + group(Cd-Cd(CCO)H) + radical(AllylJ2_triplet) + radical(Cds_S)"""),
)

species(
    label = 'C=[C]C=C1[CH]O[CH]1(24962)',
    structure = SMILES('C=[C]C=C1[CH]O[CH]1'),
    E0 = (482.206,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
    molecularWeight = (93.1033,'amu'),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[2.00529,0.0352864,4.26407e-06,-2.45842e-08,1.01482e-11,58075.1,22.2412], Tmin=(100,'K'), Tmax=(1059.43,'K')), NASAPolynomial(coeffs=[9.19063,0.0271823,-1.1198e-05,2.09564e-09,-1.47407e-13,55485,-17.8776], Tmin=(1059.43,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(482.206,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(278.535,'J/(mol*K)'), comment="""Thermo group additivity estimation: group(O2s-CsCs) + group(Cs-(Cds-Cds)OsHH) + group(Cs-(Cds-Cds)OsHH) + group(Cds-CdsCsCs) + group(Cds-Cds(Cds-Cds)H) + group(Cds-Cds(Cds-Cds)H) + group(Cds-CdsHH) + ring(Cyclobutane) + radical(C=CCJ(O)C) + radical(C=CCJ(O)C) + radical(C=CJC=C)"""),
)

species(
    label = 'C=C1[CH]C(C=O)=C1(24894)',
    structure = SMILES('C=C1[CH]C(C=O)=C1'),
    E0 = (288.82,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
    molecularWeight = (93.1033,'amu'),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[1.74733,0.0328327,2.8736e-05,-6.42855e-08,2.72716e-11,34833.1,21.0277], Tmin=(100,'K'), Tmax=(983.001,'K')), NASAPolynomial(coeffs=[16.2296,0.0146261,-5.62466e-06,1.16252e-09,-9.17833e-14,30018.3,-58.6005], Tmin=(983.001,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(288.82,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(278.535,'J/(mol*K)'), comment="""Thermo group additivity estimation: group(Cs-(Cds-Cds)(Cds-Cds)HH) + group(Cds-Cds(Cds-Cds)Cs) + group(Cd-CdCs(CO)) + group(Cds-Cds(Cds-Cds)H) + group(Cds-CdsHH) + group(Cds-O2d(Cds-Cds)H) + ring(Cyclobutene) + radical(CCJCC=O)"""),
)

species(
    label = '[CH]=C=CO[CH][C]=C(21927)',
    structure = SMILES('[CH]=C=CO[CH][C]=C'),
    E0 = (567.888,'kJ/mol'),
    modes = [
        HarmonicOscillator(frequencies=([3120,650,792.5,1650,3025,407.5,1350,352.5,540,610,2055,2950,3100,1380,975,1025,1650,1685,370,3010,987.5,1337.5,450,1655,180,180,180],'cm^-1')),
        HinderedRotor(inertia=(1.45185,'amu*angstrom^2'), symmetry=1, barrier=(33.3809,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(1.45143,'amu*angstrom^2'), symmetry=1, barrier=(33.3713,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(1.453,'amu*angstrom^2'), symmetry=1, barrier=(33.4074,'kJ/mol'), semiclassical=False),
    ],
    spinMultiplicity = 4,
    opticalIsomers = 1,
    molecularWeight = (93.1033,'amu'),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[0.391471,0.0676585,-6.93871e-05,3.5265e-08,-6.84419e-12,68441.1,27.6356], Tmin=(100,'K'), Tmax=(1379.2,'K')), NASAPolynomial(coeffs=[17.8611,0.0117575,-2.89659e-06,3.73308e-10,-2.07148e-14,64120.2,-60.4568], Tmin=(1379.2,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(567.888,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(270.22,'J/(mol*K)'), comment="""Thermo group additivity estimation: group(O2s-Cs(Cds-Cd)) + group(Cs-(Cds-Cds)OsHH) + group(Cds-CdsCsH) + group(Cds-CdsOsH) + group(Cds-CdsHH) + group(Cds-CdsHH) + group(Cdd-CdsCds) + radical(C=CCJ(O)C) + radical(C=C=CJ) + radical(Cds_S)"""),
)

species(
    label = '[C]C(C=O)=C[C]=C(24963)',
    structure = SMILES('[C]C(C=O)=C[C]=C'),
    E0 = (795.563,'kJ/mol'),
    modes = [
        HarmonicOscillator(frequencies=([1685,370,2782.5,750,1395,475,1775,1000,2950,3100,1380,975,1025,1650,350,440,435,1725,3010,987.5,1337.5,450,1655,180,180],'cm^-1')),
        HinderedRotor(inertia=(1.23283,'amu*angstrom^2'), symmetry=1, barrier=(28.3451,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(1.22773,'amu*angstrom^2'), symmetry=1, barrier=(28.2279,'kJ/mol'), semiclassical=False),
    ],
    spinMultiplicity = 5,
    opticalIsomers = 1,
    molecularWeight = (92.0954,'amu'),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[1.17334,0.0653701,-7.77773e-05,4.85088e-08,-1.22064e-11,95783.3,20.6338], Tmin=(100,'K'), Tmax=(961.607,'K')), NASAPolynomial(coeffs=[11.4795,0.0224997,-1.09043e-05,2.14692e-09,-1.53226e-13,93801.2,-28.6839], Tmin=(961.607,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(795.563,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(249.434,'J/(mol*K)'), comment="""Thermo group additivity estimation: group(Cs-(Cds-Cds)HHH) + group(Cd-CdCs(CO)) + group(Cds-Cds(Cds-Cds)H) + group(Cds-Cds(Cds-Cds)H) + group(Cds-O2d(Cds-Cds)H) + group(Cds-CdsHH) + radical(C=CJC=C) + radical(CJ3)"""),
)

species(
    label = '[CH]C(C#C[CH2])C=O(24964)',
    structure = SMILES('[CH]C(C#C[CH2])C=O'),
    E0 = (598.394,'kJ/mol'),
    modes = [
        HarmonicOscillator(frequencies=([2782.5,750,1395,475,1775,1000,1380,1390,370,380,2900,435,2100,2250,500,550,3000,3100,440,815,1455,1000,180,668.953,967.807,3967.98],'cm^-1')),
        HinderedRotor(inertia=(0.00590517,'amu*angstrom^2'), symmetry=1, barrier=(65.9828,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(5.20297,'amu*angstrom^2'), symmetry=1, barrier=(119.627,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(2.86981,'amu*angstrom^2'), symmetry=1, barrier=(65.9826,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.590213,'amu*angstrom^2'), symmetry=1, barrier=(65.9843,'kJ/mol'), semiclassical=False),
    ],
    spinMultiplicity = 4,
    opticalIsomers = 1,
    molecularWeight = (93.1033,'amu'),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[1.35491,0.0605781,-6.92373e-05,4.60238e-08,-1.26494e-11,72063.3,26.1478], Tmin=(100,'K'), Tmax=(878.415,'K')), NASAPolynomial(coeffs=[8.80035,0.0266725,-1.13367e-05,2.07856e-09,-1.41832e-13,70755.4,-8.80664], Tmin=(878.415,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(598.394,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(266.063,'J/(mol*K)'), comment="""Thermo group additivity estimation: group(Cs-(Cds-O2d)CsCsH) + group(Cs-CsHHH) + group(Cs-CtHHH) + group(Cds-OdCsH) + group(Ct-CtCs) + group(Ct-CtCs) + radical(CCJ2_triplet) + radical(Propargyl)"""),
)

species(
    label = '[CH]C([C]=O)C=C=C(24965)',
    structure = SMILES('[CH]C([C]=O)C=C=C'),
    E0 = (617.486,'kJ/mol'),
    modes = [
        HarmonicOscillator(frequencies=([1855,455,950,2950,3100,1380,975,1025,1650,540,610,2055,1380,1390,370,380,2900,435,3010,987.5,1337.5,450,1655,194.693,194.707,194.713,988.821],'cm^-1')),
        HinderedRotor(inertia=(0.0044501,'amu*angstrom^2'), symmetry=1, barrier=(0.119627,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.0044481,'amu*angstrom^2'), symmetry=1, barrier=(0.119627,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.579067,'amu*angstrom^2'), symmetry=1, barrier=(15.5789,'kJ/mol'), semiclassical=False),
    ],
    spinMultiplicity = 4,
    opticalIsomers = 1,
    molecularWeight = (93.1033,'amu'),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[1.33134,0.0622223,-7.08269e-05,4.53148e-08,-1.20017e-11,74359.5,26.1601], Tmin=(100,'K'), Tmax=(906.78,'K')), NASAPolynomial(coeffs=[9.33222,0.0269292,-1.24457e-05,2.39333e-09,-1.68362e-13,72908.4,-11.6566], Tmin=(906.78,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(617.486,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(270.22,'J/(mol*K)'), comment="""Thermo group additivity estimation: group(Cs-(Cds-O2d)CsCsH) + group(Cs-CsHHH) + group(Cds-CdsCsH) + group(Cds-OdCsH) + group(Cds-CdsHH) + group(Cdd-CdsCds) + radical(CCJ2_triplet) + radical(CC(C)CJ=O)"""),
)

species(
    label = '[CH]C(C=O)C=C=[CH](24966)',
    structure = SMILES('[CH]C(C=O)C=C=[CH]'),
    E0 = (613.274,'kJ/mol'),
    modes = [
        HarmonicOscillator(frequencies=([3120,650,792.5,1650,2782.5,750,1395,475,1775,1000,540,610,2055,1380,1390,370,380,2900,435,3010,987.5,1337.5,450,1655,300.28,303.232,303.471],'cm^-1')),
        HinderedRotor(inertia=(0.131266,'amu*angstrom^2'), symmetry=1, barrier=(8.51077,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.128883,'amu*angstrom^2'), symmetry=1, barrier=(8.48649,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.0240629,'amu*angstrom^2'), symmetry=1, barrier=(67.2931,'kJ/mol'), semiclassical=False),
    ],
    spinMultiplicity = 4,
    opticalIsomers = 1,
    molecularWeight = (93.1033,'amu'),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[1.09749,0.0668839,-8.54666e-05,6.19582e-08,-1.8247e-11,73861.5,26.4622], Tmin=(100,'K'), Tmax=(827.258,'K')), NASAPolynomial(coeffs=[9.58045,0.0258687,-1.11008e-05,2.03162e-09,-1.37866e-13,72457.9,-12.8549], Tmin=(827.258,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(613.274,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(270.22,'J/(mol*K)'), comment="""Thermo group additivity estimation: group(Cs-(Cds-O2d)CsCsH) + group(Cs-CsHHH) + group(Cds-CdsCsH) + group(Cds-OdCsH) + group(Cds-CdsHH) + group(Cdd-CdsCds) + radical(CCJ2_triplet) + radical(C=C=CJ)"""),
)

species(
    label = '[CH]C1(C=O)C=[C]C1(24967)',
    structure = SMILES('[CH]C1(C=O)C=[C]C1'),
    E0 = (695.742,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
    molecularWeight = (93.1033,'amu'),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[1.40235,0.0492961,-2.97888e-05,2.86615e-09,2.45997e-12,83778.7,23.3532], Tmin=(100,'K'), Tmax=(1063.57,'K')), NASAPolynomial(coeffs=[12.9162,0.0199037,-7.95334e-06,1.47612e-09,-1.03739e-13,80542.8,-36.6025], Tmin=(1063.57,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(695.742,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(274.378,'J/(mol*K)'), comment="""Thermo group additivity estimation: group(Cs-(Cds-O2d)CsCsCs) + group(Cs-(Cds-Cds)CsHH) + group(Cs-CsHHH) + group(Cds-CdsCsH) + group(Cds-CdsCsH) + group(Cds-OdCsH) + ring(Cyclobutene) + radical(CCJ2_triplet) + radical(cyclobutene-vinyl)"""),
)

species(
    label = '[CH]=C1C([O])C1[C]=C(24968)',
    structure = SMILES('[CH]=C1C([O])C1[C]=C'),
    E0 = (801.138,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
    molecularWeight = (93.1033,'amu'),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[0.947205,0.0595401,-5.48209e-05,2.53532e-08,-4.61065e-12,96470.9,24.0367], Tmin=(100,'K'), Tmax=(1336.35,'K')), NASAPolynomial(coeffs=[15.4105,0.0162479,-6.22675e-06,1.11075e-09,-7.54215e-14,92605.4,-49.9338], Tmin=(1336.35,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(801.138,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(278.535,'J/(mol*K)'), comment="""Thermo group additivity estimation: group(O2s-CsH) + group(Cs-(Cds-Cds)(Cds-Cds)CsH) + group(Cs-(Cds-Cds)CsOsH) + group(Cds-CdsCsCs) + group(Cds-CdsCsH) + group(Cds-CdsHH) + group(Cds-CdsHH) + ring(Methylene_cyclopropane) + radical(CC(C)OJ) + radical(Cds_P) + radical(Cds_S)"""),
)

species(
    label = 'C=[C]C1[CH]C1=C[O](24969)',
    structure = SMILES('C=[C]C1[CH]C1=C[O]'),
    E0 = (570.271,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
    molecularWeight = (93.1033,'amu'),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[1.2576,0.0475757,-8.87191e-06,-3.01789e-08,1.70641e-11,68698,20.2945], Tmin=(100,'K'), Tmax=(952.449,'K')), NASAPolynomial(coeffs=[16.8363,0.0132954,-3.93465e-06,6.98091e-10,-5.22848e-14,64317.8,-61.5208], Tmin=(952.449,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(570.271,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(278.535,'J/(mol*K)'), comment="""Thermo group additivity estimation: group(O2s-(Cds-Cd)H) + group(Cs-(Cds-Cds)(Cds-Cds)CsH) + group(Cs-(Cds-Cds)CsHH) + group(Cds-CdsCsCs) + group(Cds-CdsCsH) + group(Cds-CdsOsH) + group(Cds-CdsHH) + ring(Methylene_cyclopropane) + radical(C=COJ) + radical(Allyl_S) + radical(Cds_S)"""),
)

species(
    label = '[CH][O](751)',
    structure = SMILES('[CH][O]'),
    E0 = (424.848,'kJ/mol'),
    modes = [
        HarmonicOscillator(frequencies=([815.726,815.726,3402.81],'cm^-1')),
    ],
    spinMultiplicity = 4,
    opticalIsomers = 1,
    molecularWeight = (29.018,'amu'),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[3.86392,-0.000399472,1.49306e-05,-2.12194e-08,8.78636e-12,51105.5,7.21901], Tmin=(100,'K'), Tmax=(905.857,'K')), NASAPolynomial(coeffs=[5.97079,-0.000856178,1.03779e-06,-2.14004e-10,1.3909e-14,50360.8,-4.74054], Tmin=(905.857,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(424.848,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(58.2013,'J/(mol*K)'), comment="""Thermo group additivity estimation: group(O2s-CsH) + group(Cs-OsHHH) + radical(H3COJ) + radical(OsCsJ2H_triplet)"""),
)

species(
    label = 'C#CC=C=C(24970)',
    structure = SMILES('C#CC=C=C'),
    E0 = (411.286,'kJ/mol'),
    modes = [
        HarmonicOscillator(frequencies=([2950,3100,1380,975,1025,1650,750,770,3400,2100,3010,987.5,1337.5,450,1655,2175,525,540,610,2055],'cm^-1')),
        HinderedRotor(inertia=(1.6008,'amu*angstrom^2'), symmetry=1, barrier=(36.8056,'kJ/mol'), semiclassical=False),
    ],
    spinMultiplicity = 1,
    opticalIsomers = 1,
    molecularWeight = (64.0853,'amu'),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[2.20027,0.0319891,-8.52722e-06,-1.42168e-08,8.26989e-12,49537.8,14.2088], Tmin=(100,'K'), Tmax=(986.619,'K')), NASAPolynomial(coeffs=[11.4365,0.0119042,-4.38633e-06,8.20522e-10,-5.97198e-14,46870.3,-34.5086], Tmin=(986.619,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(411.286,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(203.705,'J/(mol*K)'), comment="""Thermo group additivity estimation: group(Cds-CdsCtH) + group(Cds-CdsHH) + group(Cdd-CdsCds) + group(Ct-Ct(Cds-Cds)) + group(Ct-CtH)"""),
)

species(
    label = '[CH]=C(C#C[CH2])C[O](24971)',
    structure = SMILES('[CH]=C(C#C[CH2])C[O]'),
    E0 = (664.143,'kJ/mol'),
    modes = [
        HarmonicOscillator(frequencies=([3120,650,792.5,1650,2750,2850,1437.5,1250,1305,750,350,350,440,435,1725,2100,2250,500,550,3000,3100,440,815,1455,1000,315.602,316.168],'cm^-1')),
        HinderedRotor(inertia=(0.00167524,'amu*angstrom^2'), symmetry=1, barrier=(0.119627,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.00168693,'amu*angstrom^2'), symmetry=1, barrier=(0.119627,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.00170209,'amu*angstrom^2'), symmetry=1, barrier=(0.119627,'kJ/mol'), semiclassical=False),
    ],
    spinMultiplicity = 4,
    opticalIsomers = 1,
    molecularWeight = (93.1033,'amu'),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[1.63952,0.056932,-7.37104e-05,6.35491e-08,-2.31134e-11,79958.1,26.4315], Tmin=(100,'K'), Tmax=(784.931,'K')), NASAPolynomial(coeffs=[4.82122,0.033696,-1.58873e-05,3.04065e-09,-2.11443e-13,79674.9,13.2301], Tmin=(784.931,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(664.143,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(270.22,'J/(mol*K)'), comment="""Thermo group additivity estimation: group(O2s-CsH) + group(Cs-(Cds-Cds)OsHH) + group(Cs-CtHHH) + group(Cds-CdsCtCs) + group(Cds-CdsHH) + group(Ct-CtCs) + group(Ct-Ct(Cds-Cds)) + radical(Cds_P) + radical(CCOJ) + radical(Propargyl)"""),
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
    label = 'C=C1[CH]C([CH]1)=C[O](24885)',
    structure = SMILES('C=C1[CH]C([CH]1)=C[O]'),
    E0 = (332.256,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
    molecularWeight = (93.1033,'amu'),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[2.34711,-0.00162064,0.000164729,-2.32489e-07,9.49087e-11,40056,19.567], Tmin=(100,'K'), Tmax=(916.82,'K')), NASAPolynomial(coeffs=[25.3663,-0.00387703,7.79823e-06,-1.58024e-09,9.57896e-14,31709.1,-111.99], Tmin=(916.82,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(332.256,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(282.692,'J/(mol*K)'), comment="""Thermo group additivity estimation: group(O2s-(Cds-Cd)H) + group(Cs-(Cds-Cds)(Cds-Cds)HH) + group(Cs-(Cds-Cds)(Cds-Cds)HH) + group(Cds-CdsCsCs) + group(Cds-CdsCsCs) + group(Cds-CdsOsH) + group(Cds-CdsHH) + ring(Cyclobutane) + radical(C=CCJC=C) + radical(C=COJ) + radical(C=CCJC=C)"""),
)

species(
    label = '[O]C=C1[CH]C[C]=C1(24798)',
    structure = SMILES('[O]C=C1[CH]C[C]=C1'),
    E0 = (428.246,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
    molecularWeight = (93.1033,'amu'),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[1.87203,0.0252799,6.10855e-05,-1.05735e-07,4.448e-11,51602.4,19.5098], Tmin=(100,'K'), Tmax=(941.537,'K')), NASAPolynomial(coeffs=[18.2956,0.0101021,-1.71259e-06,3.15951e-10,-3.16067e-14,46089.8,-71.5861], Tmin=(941.537,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(428.246,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(282.692,'J/(mol*K)'), comment="""Thermo group additivity estimation: group(O2s-(Cds-Cd)H) + group(Cs-(Cds-Cds)CsHH) + group(Cs-(Cds-Cds)CsHH) + group(Cds-Cds(Cds-Cds)Cs) + group(Cds-CdsCsH) + group(Cds-Cds(Cds-Cds)H) + group(Cds-CdsOsH) + ring(3-Methylenecyclopentene) + radical(cyclopentene-vinyl) + radical(Allyl_S) + radical(C=COJ)"""),
)

species(
    label = 'C#CC([O])[CH][C]=C(22613)',
    structure = SMILES('C#CC([O])[CH][C]=C'),
    E0 = (641.963,'kJ/mol'),
    modes = [
        HarmonicOscillator(frequencies=([750,770,3400,2100,3025,407.5,1350,352.5,2950,3100,1380,975,1025,1650,1380,1390,370,380,2900,435,1685,370,2175,525,180,180,180],'cm^-1')),
        HinderedRotor(inertia=(1.75942,'amu*angstrom^2'), symmetry=1, barrier=(40.4526,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(1.75864,'amu*angstrom^2'), symmetry=1, barrier=(40.4346,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(1.75874,'amu*angstrom^2'), symmetry=1, barrier=(40.4369,'kJ/mol'), semiclassical=False),
    ],
    spinMultiplicity = 4,
    opticalIsomers = 1,
    molecularWeight = (93.1033,'amu'),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[0.689894,0.0660626,-6.29884e-05,2.74046e-08,-3.732e-12,77335.4,25.2859], Tmin=(100,'K'), Tmax=(1001.2,'K')), NASAPolynomial(coeffs=[16.1565,0.0160403,-5.67878e-06,9.85884e-10,-6.72057e-14,73648.4,-52.2963], Tmin=(1001.2,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(641.963,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(270.22,'J/(mol*K)'), comment="""Thermo group additivity estimation: group(O2s-CsH) + group(Cs-(Cds-Cds)CsHH) + group(Cs-(Cds-Cds)CsOsH) + group(Cds-CdsCsH) + group(Cds-CdsHH) + group(Ct-CtCs) + group(Ct-CtH) + radical(CC(C)OJ) + radical(Cds_S) + radical(C=CCJCO)"""),
)

species(
    label = 'C=C=CC1=CC1[O](24972)',
    structure = SMILES('C=C=CC1=CC1[O]'),
    E0 = (517.964,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
    molecularWeight = (93.1033,'amu'),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[0.998091,0.0539889,-2.49995e-05,-1.47695e-08,1.16722e-11,62415.7,22.011], Tmin=(100,'K'), Tmax=(969.027,'K')), NASAPolynomial(coeffs=[17.8966,0.012346,-4.05392e-06,7.57962e-10,-5.73389e-14,57820.8,-65.7933], Tmin=(969.027,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(517.964,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(278.535,'J/(mol*K)'), comment="""Thermo group additivity estimation: group(O2s-CsH) + group(Cs-(Cds-Cds)(Cds-Cds)OsH) + group(Cds-Cds(Cds-Cds)Cs) + group(Cds-CdsCsH) + group(Cds-Cds(Cds-Cds)H) + group(Cds-CdsHH) + group(Cdd-CdsCds) + ring(Cyclopropene) + radical(CC(C)OJ)"""),
)

species(
    label = '[CH]=C(C=O)C=C=C(22622)',
    structure = SMILES('[CH]=C(C=O)C=C=C'),
    E0 = (367.643,'kJ/mol'),
    modes = [
        HarmonicOscillator(frequencies=([3120,650,792.5,1650,350,440,435,1725,540,610,2055,2950,3100,1380,975,1025,1650,2782.5,750,1395,475,1775,1000,3010,987.5,1337.5,450,1655],'cm^-1')),
        HinderedRotor(inertia=(1.1624,'amu*angstrom^2'), symmetry=1, barrier=(26.7259,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(1.16371,'amu*angstrom^2'), symmetry=1, barrier=(26.7559,'kJ/mol'), semiclassical=False),
    ],
    spinMultiplicity = 2,
    opticalIsomers = 1,
    molecularWeight = (93.1033,'amu'),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[0.615472,0.0660007,-6.44341e-05,3.02928e-08,-5.52751e-12,44346.1,22.6515], Tmin=(100,'K'), Tmax=(1339.17,'K')), NASAPolynomial(coeffs=[18.1439,0.0136447,-5.79027e-06,1.09868e-09,-7.74718e-14,39651.4,-67.032], Tmin=(1339.17,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(367.643,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(274.378,'J/(mol*K)'), comment="""Thermo group additivity estimation: group(Cds-Cds(Cds-O2d)(Cds-Cds)) + group(Cds-Cds(Cds-Cds)H) + group(Cds-O2d(Cds-Cds)H) + group(Cds-CdsHH) + group(Cds-CdsHH) + group(Cdd-CdsCds) + radical(Cds_P)"""),
)

species(
    label = '[CH]C(=C=O)C[C]=C(24291)',
    structure = SMILES('[CH]C(=C=O)C[C]=C'),
    E0 = (580.356,'kJ/mol'),
    modes = [
        HarmonicOscillator(frequencies=([1685,370,2750,2850,1437.5,1250,1305,750,350,2120,512.5,787.5,350,440,435,1725,2950,3100,1380,975,1025,1650,373.052,373.116,373.555,373.689,374.161],'cm^-1')),
        HinderedRotor(inertia=(0.532085,'amu*angstrom^2'), symmetry=1, barrier=(52.7173,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.530318,'amu*angstrom^2'), symmetry=1, barrier=(52.7176,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.53127,'amu*angstrom^2'), symmetry=1, barrier=(52.7155,'kJ/mol'), semiclassical=False),
    ],
    spinMultiplicity = 4,
    opticalIsomers = 1,
    molecularWeight = (93.1033,'amu'),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[1.03706,0.0712947,-9.76413e-05,8.32972e-08,-2.90902e-11,69901.5,25.3586], Tmin=(100,'K'), Tmax=(822.106,'K')), NASAPolynomial(coeffs=[5.95871,0.0372947,-1.72623e-05,3.2408e-09,-2.21901e-13,69432,4.64476], Tmin=(822.106,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(580.356,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(270.22,'J/(mol*K)'), comment="""Thermo group additivity estimation: group(Cs-Cd(CCO)HH) + group(Cs-(Cds-Cdd-O2d)HHH) + group(Cds-(Cdd-O2d)CsCs) + other(ketene_1C-C_1C-H) + group(Cds-CdsCsH) + group(Cds-CdsHH) + radical(Cds_S) + radical(AllylJ2_triplet)"""),
)

species(
    label = '[CH]=[C]CC(=[CH])C=O(22627)',
    structure = SMILES('[CH]=[C]CC(=[CH])C=O'),
    E0 = (689.978,'kJ/mol'),
    modes = [
        HarmonicOscillator(frequencies=([1685,370,2750,2850,1437.5,1250,1305,750,350,350,440,435,1725,3115,3125,620,680,785,800,1600,1700,2782.5,750,1395,475,1775,1000],'cm^-1')),
        HinderedRotor(inertia=(0.57051,'amu*angstrom^2'), symmetry=1, barrier=(13.1171,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.567077,'amu*angstrom^2'), symmetry=1, barrier=(13.0382,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.569279,'amu*angstrom^2'), symmetry=1, barrier=(13.0888,'kJ/mol'), semiclassical=False),
    ],
    spinMultiplicity = 4,
    opticalIsomers = 1,
    molecularWeight = (93.1033,'amu'),
    collisionModel = TransportData(shapeIndex=2, epsilon=(3600.17,'J/mol'), sigma=(5.90843,'angstroms'), dipoleMoment=(0,'C*m'), polarizability=(0,'angstroms^3'), rotrelaxcollnum=0, comment="""Epsilon & sigma estimated with Tc=562.34 K, Pc=39.61 bar (from Joback method)"""),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[1.42973,0.0625428,-7.42738e-05,5.1316e-08,-1.51317e-11,83072.4,25.9754], Tmin=(100,'K'), Tmax=(806.629,'K')), NASAPolynomial(coeffs=[7.73673,0.0312673,-1.61149e-05,3.24916e-09,-2.34448e-13,82055,-3.09699], Tmin=(806.629,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(689.978,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(270.22,'J/(mol*K)'), comment="""Thermo group additivity estimation: group(Cs-(Cds-Cds)(Cds-Cds)HH) + group(Cd-CdCs(CO)) + group(Cds-CdsCsH) + group(Cds-O2d(Cds-Cds)H) + group(Cds-CdsHH) + group(Cds-CdsHH) + radical(Cds_S) + radical(Cds_P) + radical(Cds_P)"""),
)

species(
    label = '[CH]=[C]C([C]=C)C=O(24973)',
    structure = SMILES('[CH]=[C]C([C]=C)C=O'),
    E0 = (695.65,'kJ/mol'),
    modes = [
        HarmonicOscillator(frequencies=([1670,1700,300,440,2782.5,750,1395,475,1775,1000,1380,1390,370,380,2900,435,3120,650,792.5,1650,2950,3100,1380,975,1025,1650,180],'cm^-1')),
        HinderedRotor(inertia=(0.394083,'amu*angstrom^2'), symmetry=1, barrier=(9.06075,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.393976,'amu*angstrom^2'), symmetry=1, barrier=(9.05829,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.393959,'amu*angstrom^2'), symmetry=1, barrier=(9.05789,'kJ/mol'), semiclassical=False),
    ],
    spinMultiplicity = 4,
    opticalIsomers = 1,
    molecularWeight = (93.1033,'amu'),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[1.11867,0.0721537,-0.000117865,1.08975e-07,-3.91257e-11,83762.7,27.2903], Tmin=(100,'K'), Tmax=(844.831,'K')), NASAPolynomial(coeffs=[5.51037,0.0330928,-1.60775e-05,3.05907e-09,-2.09531e-13,83672.6,10.7017], Tmin=(844.831,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(695.65,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(270.22,'J/(mol*K)'), comment="""Thermo group additivity estimation: group(Cs-(Cds-O2d)CsCsH) + group(Cds-CdsCsH) + group(Cds-CdsCsH) + group(Cds-OdCsH) + group(Cds-CdsHH) + group(Cds-CdsHH) + radical(Cds_P) + radical(Cds_S) + radical(Cds_S)"""),
)

species(
    label = 'C=[C]C1C=C1C=O(24974)',
    structure = SMILES('C=[C]C1C=C1C=O'),
    E0 = (457.893,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
    molecularWeight = (93.1033,'amu'),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[1.74044,0.055336,-5.67547e-05,3.52262e-08,-9.73442e-12,55148.3,20.1133], Tmin=(100,'K'), Tmax=(840.111,'K')), NASAPolynomial(coeffs=[6.47078,0.032814,-1.65429e-05,3.31687e-09,-2.39043e-13,54353.5,-1.88377], Tmin=(840.111,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(457.893,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(274.378,'J/(mol*K)'), comment="""Thermo group additivity estimation: group(Cs-(Cds-Cds)(Cds-Cds)(Cds-Cds)H) + group(Cd-CdCs(CO)) + group(Cds-CdsCsH) + group(Cds-CdsCsH) + group(Cds-O2d(Cds-Cds)H) + group(Cds-CdsHH) + ring(Cyclopropene) + radical(Cds_S)"""),
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
    E0 = (496.77,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS2',
    E0 = (692.79,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS3',
    E0 = (746.562,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS4',
    E0 = (696.285,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS5',
    E0 = (711.233,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS6',
    E0 = (712.829,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS7',
    E0 = (600.375,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS8',
    E0 = (893.246,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS9',
    E0 = (859.636,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS10',
    E0 = (872.178,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS11',
    E0 = (866.545,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS12',
    E0 = (685.636,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS13',
    E0 = (597.93,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS14',
    E0 = (853.084,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS15',
    E0 = (923.485,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS16',
    E0 = (762.063,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS17',
    E0 = (1232.77,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS18',
    E0 = (1127.67,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS19',
    E0 = (692.79,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS20',
    E0 = (726.06,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS21',
    E0 = (761.324,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS22',
    E0 = (609.974,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS23',
    E0 = (664.799,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS24',
    E0 = (875.801,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS25',
    E0 = (723.931,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS26',
    E0 = (638.097,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS27',
    E0 = (948.154,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS28',
    E0 = (663.572,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS29',
    E0 = (722.197,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS30',
    E0 = (908.882,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS31',
    E0 = (706.433,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS32',
    E0 = (915.692,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS33',
    E0 = (823.397,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS34',
    E0 = (930.83,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS35',
    E0 = (857.603,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS36',
    E0 = (1221.57,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS37',
    E0 = (612.183,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS38',
    E0 = (659.114,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS39',
    E0 = (728.021,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS40',
    E0 = (650.41,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS41',
    E0 = (621.04,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS42',
    E0 = (688.899,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS43',
    E0 = (688.899,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS44',
    E0 = (774.237,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS45',
    E0 = (614.142,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS46',
    E0 = (884.525,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS47',
    E0 = (634.825,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS48',
    E0 = (505.054,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS49',
    E0 = (881.688,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS50',
    E0 = (1007.37,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS51',
    E0 = (776.08,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS52',
    E0 = (775.685,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS53',
    E0 = (928.436,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS54',
    E0 = (695.742,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS55',
    E0 = (801.138,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS56',
    E0 = (582.023,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS57',
    E0 = (859.024,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS58',
    E0 = (806.12,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS59',
    E0 = (608.894,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS60',
    E0 = (634.825,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS61',
    E0 = (560.994,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS62',
    E0 = (805.335,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS63',
    E0 = (517.964,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS64',
    E0 = (496.77,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS65',
    E0 = (1047.14,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS66',
    E0 = (783.532,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS67',
    E0 = (835.163,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS68',
    E0 = (818.4,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS69',
    E0 = (499.699,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

reaction(
    label = 'reaction1',
    reactants = ['[CH]C(C=O)=C[C]=C(22615)'],
    products = ['C#CC=O(21959)', 'C#C[CH2](17441)'],
    transitionState = 'TS1',
    kinetics = Arrhenius(A=(5e+12,'s^-1'), n=0, Ea=(0,'kJ/mol'), T0=(1,'K'), Tmin=(300,'K'), Tmax=(1500,'K'), comment="""Exact match found for rate rule [RJJ]
Euclidian distance = 0
family: 1,4_Linear_birad_scission"""),
)

reaction(
    label = 'reaction2',
    reactants = ['[CH]C(C=O)=C[C]=C(22615)'],
    products = ['[CH]=C1[CH]C(=C)C1[O](24927)'],
    transitionState = 'TS2',
    kinetics = Arrhenius(A=(1e+10,'s^-1'), n=0, Ea=(196.02,'kJ/mol'), T0=(1,'K'), comment="""Estimated using template [R5_SD;multiplebond_intra;radadd_intra] for rate rule [R5_SD_CO;carbonylbond_intra_H;radadd_intra_cddouble]
Euclidian distance = 2.44948974278
family: Intra_R_Add_Exocyclic"""),
)

reaction(
    label = 'reaction3',
    reactants = ['H(8)', '[CH]=C(C#C[CH2])C=O(24935)'],
    products = ['[CH]C(C=O)=C[C]=C(22615)'],
    transitionState = 'TS3',
    kinetics = Arrhenius(A=(6.23e+08,'cm^3/(mol*s)'), n=1.429, Ea=(16.6816,'kJ/mol'), T0=(1,'K'), comment="""From training reaction 193 used for Ca_Ca;HJ
Exact match found for rate rule [Ca_Ca;HJ]
Euclidian distance = 0
family: R_Addition_MultipleBond"""),
)

reaction(
    label = 'reaction4',
    reactants = ['H(8)', '[CH]C(C=O)=CC#C(22634)'],
    products = ['[CH]C(C=O)=C[C]=C(22615)'],
    transitionState = 'TS4',
    kinetics = Arrhenius(A=(2.156e+09,'cm^3/(mol*s)'), n=1.502, Ea=(9.92026,'kJ/mol'), T0=(1,'K'), comment="""From training reaction 194 used for Ct-H_Ct-Cd;HJ
Exact match found for rate rule [Ct-H_Ct-Cd;HJ]
Euclidian distance = 0
family: R_Addition_MultipleBond"""),
)

reaction(
    label = 'reaction5',
    reactants = ['[CH]C(=[C]C=C)C=O(24936)'],
    products = ['[CH]C(C=O)=C[C]=C(22615)'],
    transitionState = 'TS5',
    kinetics = Arrhenius(A=(1.806e+09,'s^-1'), n=1.172, Ea=(214.463,'kJ/mol'), T0=(1,'K'), comment="""Estimated using an average for rate rule [R2H_S;Cd_rad_out_Cd;Cd_H_out_doubleC]
Euclidian distance = 0
family: intra_H_migration"""),
)

reaction(
    label = 'reaction6',
    reactants = ['[CH]C(C=O)=CC=[CH](24771)'],
    products = ['[CH]C(C=O)=C[C]=C(22615)'],
    transitionState = 'TS6',
    kinetics = Arrhenius(A=(481900,'s^-1'), n=2.375, Ea=(167.958,'kJ/mol'), T0=(1,'K'), comment="""From training reaction 121 used for R2H_D;Cd_rad_out_singleH;Cd_H_out_singleDe
Exact match found for rate rule [R2H_D;Cd_rad_out_singleH;Cd_H_out_singleDe]
Euclidian distance = 0
family: intra_H_migration"""),
)

reaction(
    label = 'reaction7',
    reactants = ['[CH]C(C=O)=C[C]=C(22615)'],
    products = ['[CH]C(=C=O)[CH]C=C(24937)'],
    transitionState = 'TS7',
    kinetics = Arrhenius(A=(1.11e+08,'s^-1'), n=1.1915, Ea=(103.605,'kJ/mol'), T0=(1,'K'), Tmin=(300,'K'), Tmax=(1600,'K'), comment="""Estimated using template [R4H_SDS;Cd_rad_out_Cd;XH_out] for rate rule [R4H_SDS;Cd_rad_out_Cd;CO_H_out]
Euclidian distance = 1.0
family: intra_H_migration"""),
)

reaction(
    label = 'reaction8',
    reactants = ['[CH]=O(373)', '[CH]=[C]C=[C][CH2](19274)'],
    products = ['[CH]C(C=O)=C[C]=C(22615)'],
    transitionState = 'TS8',
    kinetics = Arrhenius(A=(3.7839e+06,'m^3/(mol*s)'), n=0, Ea=(0,'kJ/mol'), T0=(1,'K'), comment="""Estimated using an average for rate rule [Y_rad;CO_pri_rad]
Euclidian distance = 0
family: R_Recombination
Ea raised from -3.8 to 0 kJ/mol."""),
)

reaction(
    label = 'reaction9',
    reactants = ['H(8)', '[CH]C(C#C[CH2])=C[O](24938)'],
    products = ['[CH]C(C=O)=C[C]=C(22615)'],
    transitionState = 'TS9',
    kinetics = Arrhenius(A=(4.34078e+06,'m^3/(mol*s)'), n=0.278577, Ea=(0,'kJ/mol'), T0=(1,'K'), comment="""Estimated using an average for rate rule [Y_rad;H_rad]
Euclidian distance = 0
family: R_Recombination
Ea raised from -1.4 to 0 kJ/mol."""),
)

reaction(
    label = 'reaction10',
    reactants = ['H(8)', '[CH]C([C]=O)=C[C]=C(24939)'],
    products = ['[CH]C(C=O)=C[C]=C(22615)'],
    transitionState = 'TS10',
    kinetics = Arrhenius(A=(1.53107e+07,'m^3/(mol*s)'), n=-0.133333, Ea=(0,'kJ/mol'), T0=(1,'K'), comment="""Estimated using template [CO_sec_rad;Y_rad] for rate rule [CO_rad/OneDe;H_rad]
Euclidian distance = 1.41421356237
family: R_Recombination"""),
)

reaction(
    label = 'reaction11',
    reactants = ['H(8)', '[CH]C(C=C=[CH])=C[O](22639)'],
    products = ['[CH]C(C=O)=C[C]=C(22615)'],
    transitionState = 'TS11',
    kinetics = Arrhenius(A=(4.34078e+06,'m^3/(mol*s)'), n=0.278577, Ea=(0,'kJ/mol'), T0=(1,'K'), comment="""Estimated using an average for rate rule [Y_rad;H_rad]
Euclidian distance = 0
family: R_Recombination
Ea raised from -1.4 to 0 kJ/mol."""),
)

reaction(
    label = 'reaction12',
    reactants = ['[CH]C(C=O)=C[C]=C(22615)'],
    products = ['[CH]C1([CH]C1=C)C=O(24940)'],
    transitionState = 'TS12',
    kinetics = Arrhenius(A=(5.605e+12,'s^-1'), n=0.275, Ea=(188.866,'kJ/mol'), T0=(1,'K'), comment="""Estimated using an average for rate rule [R3_D;doublebond_intra_pri;radadd_intra_cddouble]
Euclidian distance = 0
family: Intra_R_Add_Endocyclic
Ea raised from 186.0 to 188.9 kJ/mol to match endothermicity of reaction."""),
)

reaction(
    label = 'reaction13',
    reactants = ['[CH]C(C=O)=C[C]=C(22615)'],
    products = ['[CH]C1[CH]OC(=C)C=1(24910)'],
    transitionState = 'TS13',
    kinetics = Arrhenius(A=(7.96106e+11,'s^-1'), n=0.00276955, Ea=(101.161,'kJ/mol'), T0=(1,'K'), comment="""Estimated using template [R5_SD;multiplebond_intra;radadd_intra] for rate rule [R5_SD_CO;carbonyl_intra_H;radadd_intra_cddouble]
Euclidian distance = 2.44948974278
family: Intra_R_Add_Endocyclic"""),
)

reaction(
    label = 'reaction14',
    reactants = ['[CH]C([C]=O)C=[C][CH2](24941)'],
    products = ['[CH]C(C=O)=C[C]=C(22615)'],
    transitionState = 'TS14',
    kinetics = Arrhenius(A=(1.949e+11,'s^-1'), n=0.486, Ea=(22.8614,'kJ/mol'), T0=(1,'K'), comment="""Estimated using an average for rate rule [R2radExo;Y_rad;XH_Rrad]
Euclidian distance = 0
family: Intra_Disproportionation"""),
)

reaction(
    label = 'reaction15',
    reactants = ['[CH]C(=[C][C]=C)C[O](24942)'],
    products = ['[CH]C(C=O)=C[C]=C(22615)'],
    transitionState = 'TS15',
    kinetics = Arrhenius(A=(1.4874e+09,'s^-1'), n=1.045, Ea=(63.4002,'kJ/mol'), T0=(1,'K'), comment="""Estimated using an average for rate rule [R3radExo;Y_rad;XH_Rrad]
Euclidian distance = 0
Multiplied by reaction path degeneracy 2.0
family: Intra_Disproportionation"""),
)

reaction(
    label = 'reaction16',
    reactants = ['[CH]C([C]=[C][CH2])=CO(24943)'],
    products = ['[CH]C(C=O)=C[C]=C(22615)'],
    transitionState = 'TS16',
    kinetics = Arrhenius(A=(5.14222e+08,'s^-1'), n=0.311, Ea=(24.9733,'kJ/mol'), T0=(1,'K'), comment="""Estimated using template [R4;Y_rad;XH_Rrad] for rate rule [R4radEndo;Y_rad;XH_Rrad]
Euclidian distance = 1.0
family: Intra_Disproportionation"""),
)

reaction(
    label = 'reaction17',
    reactants = ['[C-]#[O+](374)', '[CH]C=C[C]=C(17762)'],
    products = ['[CH]C(C=O)=C[C]=C(22615)'],
    transitionState = 'TS17',
    kinetics = Arrhenius(A=(0.0591985,'m^3/(mol*s)'), n=2.3675, Ea=(305.306,'kJ/mol'), T0=(1,'K'), comment="""Estimated using an average for rate rule [CO;R_H]
Euclidian distance = 0
family: 1,2_Insertion_CO"""),
)

reaction(
    label = 'reaction18',
    reactants = ['[C]=C(584)', '[CH]C(=[CH])C=O(21223)'],
    products = ['[CH]C(C=O)=C[C]=C(22615)'],
    transitionState = 'TS18',
    kinetics = Arrhenius(A=(2.29708e+06,'m^3/(mol*s)'), n=0.575199, Ea=(34.3157,'kJ/mol'), T0=(1,'K'), comment="""Estimated using template [Y_rad;Birad] for rate rule [Cd_pri_rad;Birad]
Euclidian distance = 2.0
Multiplied by reaction path degeneracy 2.0
family: Birad_R_Recombination"""),
)

reaction(
    label = 'reaction19',
    reactants = ['[CH]C(C=O)=C[C]=C(22615)'],
    products = ['[CH]C1=COC1[C]=C(24944)'],
    transitionState = 'TS19',
    kinetics = Arrhenius(A=(1e+10,'s^-1'), n=0, Ea=(196.02,'kJ/mol'), T0=(1,'K'), comment="""Estimated using template [R5_SD_D;doublebond_intra;radadd_intra] for rate rule [R5_SD_D;doublebond_intra;radadd_intra_O]
Euclidian distance = 1.0
family: Intra_R_Add_Exocyclic"""),
)

reaction(
    label = 'reaction20',
    reactants = ['H(8)', '[CH]C(=C=O)C=C=C(24945)'],
    products = ['[CH]C(C=O)=C[C]=C(22615)'],
    transitionState = 'TS20',
    kinetics = Arrhenius(A=(7048.5,'m^3/(mol*s)'), n=1.43107, Ea=(2.92155,'kJ/mol'), T0=(1,'K'), Tmin=(303.03,'K'), Tmax=(2000,'K'), comment="""From training reaction 188 used for Ck_O;HJ
Exact match found for rate rule [Ck_O;HJ]
Euclidian distance = 0
family: R_Addition_MultipleBond"""),
)

reaction(
    label = 'reaction21',
    reactants = ['[CH]C(=[C]O)C=C=C(24946)'],
    products = ['[CH]C(C=O)=C[C]=C(22615)'],
    transitionState = 'TS21',
    kinetics = Arrhenius(A=(4.96519e+09,'s^-1'), n=1.05826, Ea=(162.779,'kJ/mol'), T0=(1,'K'), comment="""Estimated using template [R2H_S;Cd_rad_out_Cd;XH_out] for rate rule [R2H_S;Cd_rad_out_Cd;O_H_out]
Euclidian distance = 1.0
family: intra_H_migration"""),
)

reaction(
    label = 'reaction22',
    reactants = ['[CH]C(C#C[CH2])=CO(24947)'],
    products = ['[CH]C(C=O)=C[C]=C(22615)'],
    transitionState = 'TS22',
    kinetics = Arrhenius(A=(1.11e+08,'s^-1'), n=1.1915, Ea=(103.605,'kJ/mol'), T0=(1,'K'), comment="""Estimated using template [R4H_SDS;Cd_rad_out_double;XH_out] for rate rule [R4H_SDS;Cd_rad_out_double;O_H_out]
Euclidian distance = 1.0
family: intra_H_migration"""),
)

reaction(
    label = 'reaction23',
    reactants = ['[CH]C(C=C=[CH])=CO(24948)'],
    products = ['[CH]C(C=O)=C[C]=C(22615)'],
    transitionState = 'TS23',
    kinetics = Arrhenius(A=(1.86943e+06,'s^-1'), n=1.85754, Ea=(151.521,'kJ/mol'), T0=(1,'K'), comment="""Estimated using template [RnH;Cd_rad_out_singleH;XH_out] for rate rule [R6H;Cd_rad_out_singleH;O_H_out]
Euclidian distance = 2.2360679775
family: intra_H_migration"""),
)

reaction(
    label = 'reaction24',
    reactants = ['[CH][C]=C[O](21209)', 'C#C[CH2](17441)'],
    products = ['[CH]C(C=O)=C[C]=C(22615)'],
    transitionState = 'TS24',
    kinetics = Arrhenius(A=(40577.9,'m^3/(mol*s)'), n=0.702818, Ea=(0,'kJ/mol'), T0=(1,'K'), comment="""Estimated using an average for rate rule [Y_rad;Cd_allenic]
Euclidian distance = 0
family: R_Recombination
Ea raised from -11.4 to 0 kJ/mol."""),
)

reaction(
    label = 'reaction25',
    reactants = ['[CH]C(C=O)=C[C]=C(22615)'],
    products = ['[CH]C1([CH]O1)C=C=C(24949)'],
    transitionState = 'TS25',
    kinetics = Arrhenius(A=(2.42595e+10,'s^-1'), n=0.7335, Ea=(227.162,'kJ/mol'), T0=(1,'K'), comment="""Estimated using template [R3_D;doublebond_intra_pri;radadd_intra] for rate rule [R3_D;doublebond_intra_pri;radadd_intra_O]
Euclidian distance = 1.0
family: Intra_R_Add_Endocyclic
Ea raised from 226.1 to 227.2 kJ/mol to match endothermicity of reaction."""),
)

reaction(
    label = 'reaction26',
    reactants = ['[CH]C(C=O)=C[C]=C(22615)'],
    products = ['[CH]C1C=[C]COC=1(24822)'],
    transitionState = 'TS26',
    kinetics = Arrhenius(A=(2.69e+10,'s^-1'), n=0.239, Ea=(141.327,'kJ/mol'), T0=(1,'K'), comment="""Estimated using template [R6_SDS_D;doublebond_intra_CdCdd;radadd_intra] for rate rule [R6_SDS_D;doublebond_intra_CdCdd;radadd_intra_O]
Euclidian distance = 1.0
family: Intra_R_Add_Endocyclic"""),
)

reaction(
    label = 'reaction27',
    reactants = ['[CH]C(C#C[CH2])[CH][O](24950)'],
    products = ['[CH]C(C=O)=C[C]=C(22615)'],
    transitionState = 'TS27',
    kinetics = Arrhenius(A=(1.949e+11,'s^-1'), n=0.486, Ea=(22.8614,'kJ/mol'), T0=(1,'K'), comment="""Estimated using an average for rate rule [R2radExo;Y_rad_De;XH_Rrad]
Euclidian distance = 0
family: Intra_Disproportionation"""),
)

reaction(
    label = 'reaction28',
    reactants = ['[CH]C([C]=C[CH2])=C[O](24951)'],
    products = ['[CH]C(C=O)=C[C]=C(22615)'],
    transitionState = 'TS28',
    kinetics = Arrhenius(A=(1.949e+11,'s^-1'), n=0.486, Ea=(22.8614,'kJ/mol'), T0=(1,'K'), comment="""Estimated using an average for rate rule [R2radExo;Y_rad_De;XH_Rrad]
Euclidian distance = 0
family: Intra_Disproportionation"""),
)

reaction(
    label = 'reaction29',
    reactants = ['[CH]C=CC([CH])=C[O](24952)'],
    products = ['[CH]C(C=O)=C[C]=C(22615)'],
    transitionState = 'TS29',
    kinetics = Arrhenius(A=(1.4733e+10,'s^-1'), n=0.2847, Ea=(27.8529,'kJ/mol'), T0=(1,'K'), comment="""Estimated using average of templates [Rn;Y_rad;XH_Rrad_De] + [R2radExo;Y_rad;XH_Rrad] for rate rule [R2radExo;Y_rad;XH_Rrad_De]
Euclidian distance = 1.0
family: Intra_Disproportionation"""),
)

reaction(
    label = 'reaction30',
    reactants = ['[CH]C(=[C][O])C[C]=C(24953)'],
    products = ['[CH]C(C=O)=C[C]=C(22615)'],
    transitionState = 'TS30',
    kinetics = Arrhenius(A=(4.00798e+09,'s^-1'), n=0.37, Ea=(78.2471,'kJ/mol'), T0=(1,'K'), comment="""Estimated using average of templates [R3;Y_rad;XH_Rrad_De] + [R3radExo;Y_rad;XH_Rrad] for rate rule [R3radExo;Y_rad;XH_Rrad_De]
Euclidian distance = 1.0
Multiplied by reaction path degeneracy 2.0
family: Intra_Disproportionation"""),
)

reaction(
    label = 'reaction31',
    reactants = ['[CH]C(=[C][O])C=C[CH2](24954)'],
    products = ['[CH]C(C=O)=C[C]=C(22615)'],
    transitionState = 'TS31',
    kinetics = Arrhenius(A=(5.14222e+08,'s^-1'), n=0.311, Ea=(24.9733,'kJ/mol'), T0=(1,'K'), comment="""Estimated using template [R4;Y_rad;XH_Rrad] for rate rule [R4radEndo;Y_rad;XH_Rrad]
Euclidian distance = 1.0
family: Intra_Disproportionation"""),
)

reaction(
    label = 'reaction32',
    reactants = ['[CH][C]=CC([CH])C=O(21904)'],
    products = ['[CH]C(C=O)=C[C]=C(22615)'],
    transitionState = 'TS32',
    kinetics = Arrhenius(A=(5.14222e+08,'s^-1'), n=0.311, Ea=(24.9733,'kJ/mol'), T0=(1,'K'), comment="""Estimated using template [R4;Y_rad;XH_Rrad] for rate rule [R4radExo;Y_rad;XH_Rrad]
Euclidian distance = 1.0
family: Intra_Disproportionation"""),
)

reaction(
    label = 'reaction33',
    reactants = ['[CH][C]([C]=O)C=[C]C(24955)'],
    products = ['[CH]C(C=O)=C[C]=C(22615)'],
    transitionState = 'TS33',
    kinetics = Arrhenius(A=(6.37831e+09,'s^-1'), n=0.137, Ea=(24.9733,'kJ/mol'), T0=(1,'K'), comment="""Estimated using template [R5;Y_rad;XH_Rrad] for rate rule [R5radEndo;Y_rad;XH_Rrad]
Euclidian distance = 1.0
Multiplied by reaction path degeneracy 3.0
family: Intra_Disproportionation"""),
)

reaction(
    label = 'reaction34',
    reactants = ['[CH][C]([CH]C#C)C[O](24779)'],
    products = ['[CH]C(C=O)=C[C]=C(22615)'],
    transitionState = 'TS34',
    kinetics = Arrhenius(A=(4.25221e+09,'s^-1'), n=0.137, Ea=(24.9733,'kJ/mol'), T0=(1,'K'), comment="""Estimated using template [R5;Y_rad;XH_Rrad] for rate rule [R5radEndo;Y_rad;XH_Rrad]
Euclidian distance = 1.0
Multiplied by reaction path degeneracy 2.0
family: Intra_Disproportionation"""),
)

reaction(
    label = 'reaction35',
    reactants = ['[CH]C(=CC#C)C[O](24956)'],
    products = ['[CH]C(C=O)=C[C]=C(22615)'],
    transitionState = 'TS35',
    kinetics = Arrhenius(A=(8.3295e+09,'s^-1'), n=0.737748, Ea=(218.723,'kJ/mol'), T0=(1,'K'), comment="""Estimated using an average for rate rule [1_3_unsaturated_pentane_backbone;CH_end;CtH_2]
Euclidian distance = 0
Multiplied by reaction path degeneracy 2.0
family: Intra_ene_reaction"""),
)

reaction(
    label = 'reaction36',
    reactants = ['O(T)(63)', '[CH]C(=[CH])C=C=C(21246)'],
    products = ['[CH]C(C=O)=C[C]=C(22615)'],
    transitionState = 'TS36',
    kinetics = Arrhenius(A=(187219,'m^3/(mol*s)'), n=1.13083, Ea=(163.847,'kJ/mol'), T0=(1,'K'), Tmin=(303.03,'K'), Tmax=(2000,'K'), comment="""Estimated using template [Y_rad;O_birad] for rate rule [Cd_pri_rad;O_birad]
Euclidian distance = 2.0
Multiplied by reaction path degeneracy 2.0
family: Birad_R_Recombination"""),
)

reaction(
    label = 'reaction37',
    reactants = ['[CH]C(C=O)=C[C]=C(22615)'],
    products = ['C=[C]C=C1[CH]C1[O](24957)'],
    transitionState = 'TS37',
    kinetics = Arrhenius(A=(9.36651e+10,'s^-1'), n=0.5685, Ea=(115.413,'kJ/mol'), T0=(1,'K'), comment="""Estimated using template [R4_D;multiplebond_intra;radadd_intra_cdsingleH] for rate rule [R4_D_CO;carbonylbond_intra_H;radadd_intra_cdsingleH]
Euclidian distance = 2.2360679775
family: Intra_R_Add_Exocyclic
Ea raised from 113.5 to 115.4 kJ/mol to match endothermicity of reaction."""),
)

reaction(
    label = 'reaction38',
    reactants = ['[CH]C(C=O)=C[C]=C(22615)'],
    products = ['[CH]=C1C=[C]CC1[O](24833)'],
    transitionState = 'TS38',
    kinetics = Arrhenius(A=(2.42978e+08,'s^-1'), n=0.660014, Ea=(162.344,'kJ/mol'), T0=(1,'K'), comment="""Estimated using template [R6;multiplebond_intra;radadd_intra_cs2H] for rate rule [R6;carbonylbond_intra_H;radadd_intra_cs2H]
Euclidian distance = 2.0
family: Intra_R_Add_Exocyclic
Ea raised from 157.8 to 162.3 kJ/mol to match endothermicity of reaction."""),
)

reaction(
    label = 'reaction39',
    reactants = ['C#CC=O(21959)', '[CH][C]=C(18825)'],
    products = ['[CH]C(C=O)=C[C]=C(22615)'],
    transitionState = 'TS39',
    kinetics = Arrhenius(A=(0.0669079,'m^3/(mol*s)'), n=2.39465, Ea=(29.077,'kJ/mol'), T0=(1,'K'), comment="""Estimated using template [Ct-De_Ct-H;CJ] for rate rule [Ct-CO_Ct-H;CJ]
Euclidian distance = 1.0
family: R_Addition_MultipleBond"""),
)

reaction(
    label = 'reaction40',
    reactants = ['[CH]=O(373)', '[CH]=C=C[C]=C(19277)'],
    products = ['[CH]C(C=O)=C[C]=C(22615)'],
    transitionState = 'TS40',
    kinetics = Arrhenius(A=(0.0669079,'m^3/(mol*s)'), n=2.39465, Ea=(29.077,'kJ/mol'), T0=(1,'K'), comment="""Estimated using template [Ct-De_Ct-H;CJ] for rate rule [Ct-De_Ct-H;CO_pri_rad]
Euclidian distance = 2.0
family: R_Addition_MultipleBond"""),
)

reaction(
    label = 'reaction41',
    reactants = ['[CH]=C=C[O](8556)', 'C#C[CH2](17441)'],
    products = ['[CH]C(C=O)=C[C]=C(22615)'],
    transitionState = 'TS41',
    kinetics = Arrhenius(A=(0.523563,'m^3/(mol*s)'), n=2.10494, Ea=(22.6844,'kJ/mol'), T0=(1,'K'), comment="""Estimated using an average for rate rule [Ct_Ct;CJ]
Euclidian distance = 0
family: R_Addition_MultipleBond"""),
)

reaction(
    label = 'reaction42',
    reactants = ['[CH]C(C=O)=C[C]=C(22615)'],
    products = ['[CH2]C#CC([CH2])=C[O](24958)'],
    transitionState = 'TS42',
    kinetics = Arrhenius(A=(6718.85,'s^-1'), n=2.58467, Ea=(192.129,'kJ/mol'), T0=(1,'K'), comment="""Estimated using an average for rate rule [R3H_DS;Cd_rad_out_singleH;XH_out]
Euclidian distance = 0
family: intra_H_migration"""),
)

reaction(
    label = 'reaction43',
    reactants = ['[CH]C(C=O)=C[C]=C(22615)'],
    products = ['[CH2]C([C]=O)=C[C]=C(24959)'],
    transitionState = 'TS43',
    kinetics = Arrhenius(A=(6718.85,'s^-1'), n=2.58467, Ea=(192.129,'kJ/mol'), T0=(1,'K'), comment="""Estimated using template [R3H_DS;Cd_rad_out_singleH;XH_out] for rate rule [R3H_DS;Cd_rad_out_singleH;CO_H_out]
Euclidian distance = 1.0
family: intra_H_migration"""),
)

reaction(
    label = 'reaction44',
    reactants = ['[CH]C(C=O)=C[C]=C(22615)'],
    products = ['[CH]C(C#CC)=C[O](24960)'],
    transitionState = 'TS44',
    kinetics = Arrhenius(A=(1.59e+07,'s^-1'), n=1.4638, Ea=(277.467,'kJ/mol'), T0=(1,'K'), Tmin=(300,'K'), Tmax=(1600,'K'), comment="""Estimated using template [R3Hall;C_rad_out_2H;Cd_H_out_singleDe] for rate rule [R3HJ;C_rad_out_2H;Cd_H_out_singleDe]
Euclidian distance = 1.0
family: intra_H_migration"""),
)

reaction(
    label = 'reaction45',
    reactants = ['[CH]C(=C=O)C=[C]C(24961)'],
    products = ['[CH]C(C=O)=C[C]=C(22615)'],
    transitionState = 'TS45',
    kinetics = Arrhenius(A=(505536,'s^-1'), n=1.7378, Ea=(41.5716,'kJ/mol'), T0=(1,'K'), comment="""Estimated using template [R5Hall;Y_rad_out;Cs_H_out_2H] for rate rule [R5HJ_3;CO_rad_out;Cs_H_out_2H]
Euclidian distance = 1.41421356237
Multiplied by reaction path degeneracy 3.0
family: intra_H_migration"""),
)

reaction(
    label = 'reaction46',
    reactants = ['[CH]=C=C[O](8556)', '[CH][C]=C(18825)'],
    products = ['[CH]C(C=O)=C[C]=C(22615)'],
    transitionState = 'TS46',
    kinetics = Arrhenius(A=(1.9789e+07,'m^3/(mol*s)'), n=-0.126319, Ea=(0,'kJ/mol'), T0=(1,'K'), comment="""Estimated using an average for rate rule [Y_rad;Y_rad]
Euclidian distance = 0
family: R_Recombination
Ea raised from -15.6 to -15.6 kJ/mol.
Ea raised from -15.6 to 0 kJ/mol."""),
)

reaction(
    label = 'reaction47',
    reactants = ['[CH]C(C=O)=C[C]=C(22615)'],
    products = ['C=[C]C=C1[CH]O[CH]1(24962)'],
    transitionState = 'TS47',
    kinetics = Arrhenius(A=(1.953e+11,'s^-1'), n=0.387, Ea=(138.055,'kJ/mol'), T0=(1,'K'), comment="""Estimated using template [R4_D;multiplebond_intra;radadd_intra_cdsingleH] for rate rule [R4_D_CO;carbonyl_intra_H;radadd_intra_cdsingleH]
Euclidian distance = 2.2360679775
family: Intra_R_Add_Endocyclic"""),
)

reaction(
    label = 'reaction48',
    reactants = ['[CH]C(C=O)=C[C]=C(22615)'],
    products = ['C=C1[CH]C(C=O)=C1(24894)'],
    transitionState = 'TS48',
    kinetics = Arrhenius(A=(1.62e+12,'s^-1'), n=-0.305, Ea=(8.28432,'kJ/mol'), T0=(1,'K'), comment="""Estimated using template [R4;Y_rad_out;Ypri_rad_out] for rate rule [R4_DSD;CdsingleH_rad_out;Ypri_rad_out]
Euclidian distance = 2.2360679775
family: Birad_recombination"""),
)

reaction(
    label = 'reaction49',
    reactants = ['[CH]=C=CO[CH][C]=C(21927)'],
    products = ['[CH]C(C=O)=C[C]=C(22615)'],
    transitionState = 'TS49',
    kinetics = Arrhenius(A=(7040,'s^-1'), n=2.66, Ea=(313.8,'kJ/mol'), T0=(1,'K'), Tmin=(300,'K'), Tmax=(1500,'K'), comment="""Estimated using an average for rate rule [R_ROR;R1_doublebond;R2_doublebond_H;R_O_C]
Euclidian distance = 0
family: ketoenol"""),
)

reaction(
    label = 'reaction50',
    reactants = ['H(8)', '[C]C(C=O)=C[C]=C(24963)'],
    products = ['[CH]C(C=O)=C[C]=C(22615)'],
    transitionState = 'TS50',
    kinetics = Arrhenius(A=(1e+07,'m^3/(mol*s)'), n=0, Ea=(0,'kJ/mol'), T0=(1,'K'), comment="""Estimated using an average for rate rule [H_rad;Birad]
Euclidian distance = 0
family: Birad_R_Recombination"""),
)

reaction(
    label = 'reaction51',
    reactants = ['[CH]C(C#C[CH2])C=O(24964)'],
    products = ['[CH]C(C=O)=C[C]=C(22615)'],
    transitionState = 'TS51',
    kinetics = Arrhenius(A=(3.9024e+09,'s^-1'), n=1.09833, Ea=(177.686,'kJ/mol'), T0=(1,'K'), comment="""Estimated using an average for rate rule [R2H_S;Cd_rad_out_double;XH_out]
Euclidian distance = 0
family: intra_H_migration"""),
)

reaction(
    label = 'reaction52',
    reactants = ['[CH]C([C]=O)C=C=C(24965)'],
    products = ['[CH]C(C=O)=C[C]=C(22615)'],
    transitionState = 'TS52',
    kinetics = Arrhenius(A=(4.89315e+07,'s^-1'), n=1.5961, Ea=(158.199,'kJ/mol'), T0=(1,'K'), comment="""Estimated using template [R2H_S;Y_rad_out;XH_out] for rate rule [R2H_S;CO_rad_out;XH_out]
Euclidian distance = 1.0
family: intra_H_migration"""),
)

reaction(
    label = 'reaction53',
    reactants = ['[CH]C(C=O)C=C=[CH](24966)'],
    products = ['[CH]C(C=O)=C[C]=C(22615)'],
    transitionState = 'TS53',
    kinetics = Arrhenius(A=(3293.16,'s^-1'), n=2.5965, Ea=(315.162,'kJ/mol'), T0=(1,'K'), comment="""Estimated using an average for rate rule [R4H_MMS;Cd_rad_out_singleH;XH_out]
Euclidian distance = 0
family: intra_H_migration"""),
)

reaction(
    label = 'reaction54',
    reactants = ['[CH]C(C=O)=C[C]=C(22615)'],
    products = ['[CH]C1(C=O)C=[C]C1(24967)'],
    transitionState = 'TS54',
    kinetics = Arrhenius(A=(3.01156e+10,'s^-1'), n=0.428741, Ea=(198.973,'kJ/mol'), T0=(1,'K'), comment="""Estimated using average of templates [Rn;doublebond_intra_CdCdd;radadd_intra_cs] + [R4;doublebond_intra;radadd_intra_cs] for rate rule [R4;doublebond_intra_CdCdd;radadd_intra_cs]
Euclidian distance = 1.0
family: Intra_R_Add_Endocyclic
Ea raised from 195.9 to 199.0 kJ/mol to match endothermicity of reaction."""),
)

reaction(
    label = 'reaction55',
    reactants = ['[CH]C(C=O)=C[C]=C(22615)'],
    products = ['[CH]=C1C([O])C1[C]=C(24968)'],
    transitionState = 'TS55',
    kinetics = Arrhenius(A=(7.18499e+11,'s^-1'), n=-0.0609598, Ea=(304.368,'kJ/mol'), T0=(1,'K'), comment="""Estimated using template [R4_S_D;doublebond_intra;radadd_intra_cs] for rate rule [R4_S_(Cd)_D;doublebond_intra;radadd_intra_cs]
Euclidian distance = 2.0
family: Intra_R_Add_Exocyclic
Ea raised from 302.6 to 304.4 kJ/mol to match endothermicity of reaction."""),
)

reaction(
    label = 'reaction56',
    reactants = ['[CH]C(C=O)=C[C]=C(22615)'],
    products = ['C=[C]C1[CH]C1=C[O](24969)'],
    transitionState = 'TS56',
    kinetics = Arrhenius(A=(2.54e+10,'s^-1'), n=0.69, Ea=(85.2532,'kJ/mol'), T0=(1,'K'), comment="""Estimated using an average for rate rule [R4_D_D;doublebond_intra;radadd_intra_cdsingleH]
Euclidian distance = 0
family: Intra_R_Add_Exocyclic"""),
)

reaction(
    label = 'reaction57',
    reactants = ['[CH][O](751)', 'C#CC=C=C(24970)'],
    products = ['[CH]C(C=O)=C[C]=C(22615)'],
    transitionState = 'TS57',
    kinetics = Arrhenius(A=(3.01923,'m^3/(mol*s)'), n=1.94267, Ea=(22.8894,'kJ/mol'), T0=(1,'K'), comment="""Estimated using template [Ct-Cd_Ct-H;YJ] for rate rule [Ct-Cd_Ct-H;Y_1centerbirad]
Euclidian distance = 1.0
family: R_Addition_MultipleBond"""),
)

reaction(
    label = 'reaction58',
    reactants = ['[CH]=C(C#C[CH2])C[O](24971)'],
    products = ['[CH]C(C=O)=C[C]=C(22615)'],
    transitionState = 'TS58',
    kinetics = Arrhenius(A=(4.823e+09,'s^-1'), n=1.00333, Ea=(141.977,'kJ/mol'), T0=(1,'K'), comment="""Estimated using template [R3H_SS;Cd_rad_out_double;XH_out] for rate rule [R3H_SS_2Cd;Cd_rad_out_double;XH_out]
Euclidian distance = 1.0
Multiplied by reaction path degeneracy 2.0
family: intra_H_migration"""),
)

reaction(
    label = 'reaction59',
    reactants = ['[CH]C(C=O)=C[C]=C(22615)'],
    products = ['[CH]=C=CC([CH2])=C[O](21884)'],
    transitionState = 'TS59',
    kinetics = Arrhenius(A=(1.99006e+08,'s^-1'), n=1.38995, Ea=(112.125,'kJ/mol'), T0=(1,'K'), comment="""Estimated using average of templates [RnH;Cd_rad_out_singleH;Cd_H_out_singleH] + [R5H;Cd_rad_out_singleH;XH_out] for rate rule [R5H;Cd_rad_out_singleH;Cd_H_out_singleH]
Euclidian distance = 2.0
Multiplied by reaction path degeneracy 2.0
family: intra_H_migration"""),
)

reaction(
    label = 'reaction60',
    reactants = ['[CH]C(C=O)=C[C]=C(22615)'],
    products = ['C=C1[CH]C([CH]1)=C[O](24885)'],
    transitionState = 'TS60',
    kinetics = Arrhenius(A=(1.953e+11,'s^-1'), n=0.387, Ea=(138.055,'kJ/mol'), T0=(1,'K'), comment="""Estimated using an average for rate rule [R4_D_D;doublebond_intra;radadd_intra_cdsingleH]
Euclidian distance = 0
family: Intra_R_Add_Endocyclic"""),
)

reaction(
    label = 'reaction61',
    reactants = ['[CH]C(C=O)=C[C]=C(22615)'],
    products = ['[O]C=C1[CH]C[C]=C1(24798)'],
    transitionState = 'TS61',
    kinetics = Arrhenius(A=(7.58e+12,'s^-1'), n=-0.292, Ea=(64.2244,'kJ/mol'), T0=(1,'K'), comment="""From training reaction 838 used for R5_MS;doublebond_intra_CdCdd;radadd_intra_cdsingleH
Exact match found for rate rule [R5_MS;doublebond_intra_CdCdd;radadd_intra_cdsingleH]
Euclidian distance = 0
family: Intra_R_Add_Endocyclic"""),
)

reaction(
    label = 'reaction62',
    reactants = ['C#CC([O])[CH][C]=C(22613)'],
    products = ['[CH]C(C=O)=C[C]=C(22615)'],
    transitionState = 'TS62',
    kinetics = Arrhenius(A=(3.36138e+09,'s^-1'), n=0.97875, Ea=(163.373,'kJ/mol'), T0=(1,'K'), comment="""Estimated using an average for rate rule [cCsCJ;CJ;C]
Euclidian distance = 0
family: 1,2_shiftC"""),
)

reaction(
    label = 'reaction63',
    reactants = ['[CH]C(C=O)=C[C]=C(22615)'],
    products = ['C=C=CC1=CC1[O](24972)'],
    transitionState = 'TS63',
    kinetics = Arrhenius(A=(2.18e+16,'s^-1'), n=0, Ea=(21.1947,'kJ/mol'), T0=(1,'K'), comment="""Estimated using template [R3;Y_rad_out;Ypri_rad_out] for rate rule [R3_SD;Y_rad_out;CdsinglepriH_rad_out]
Euclidian distance = 2.2360679775
family: Birad_recombination
Ea raised from 19.1 to 21.2 kJ/mol to match endothermicity of reaction."""),
)

reaction(
    label = 'reaction64',
    reactants = ['[CH]C(C=O)=C[C]=C(22615)'],
    products = ['[CH]=C(C=O)C=C=C(22622)'],
    transitionState = 'TS64',
    kinetics = Arrhenius(A=(1e+10,'s^-1'), n=0, Ea=(0,'kJ/mol'), T0=(1,'K'), comment="""Estimated using an average for rate rule [NOS]
Euclidian distance = 0
family: 1,2-Birad_to_alkene"""),
)

reaction(
    label = 'reaction65',
    reactants = ['[CH][O](751)', '[CH]=C=C[C]=C(19277)'],
    products = ['[CH]C(C=O)=C[C]=C(22615)'],
    transitionState = 'TS65',
    kinetics = Arrhenius(A=(1.14854e+06,'m^3/(mol*s)'), n=0.575199, Ea=(34.3157,'kJ/mol'), T0=(1,'K'), comment="""Estimated using an average for rate rule [Y_rad;Birad]
Euclidian distance = 0
family: Birad_R_Recombination"""),
)

reaction(
    label = 'reaction66',
    reactants = ['[CH]C(=C=O)C[C]=C(24291)'],
    products = ['[CH]C(C=O)=C[C]=C(22615)'],
    transitionState = 'TS66',
    kinetics = Arrhenius(A=(9.59054e+08,'s^-1'), n=1.35081, Ea=(203.176,'kJ/mol'), T0=(1,'K'), comment="""Estimated using template [R3H_SS_2Cd;Y_rad_out;XH_out] for rate rule [R3H_SS_2Cd;CO_rad_out;XH_out]
Euclidian distance = 1.0
Multiplied by reaction path degeneracy 2.0
family: intra_H_migration"""),
)

reaction(
    label = 'reaction12',
    reactants = ['[CH]=[C]CC(=[CH])C=O(22627)'],
    products = ['[CH]C(C=O)=C[C]=C(22615)'],
    transitionState = 'TS67',
    kinetics = Arrhenius(A=(1.846e+10,'s^-1'), n=0.74, Ea=(145.185,'kJ/mol'), T0=(1,'K'), comment="""Estimated using template [R3Hall;Cd_rad_out_singleH;Cs_H_out_1H] for rate rule [R3HJ;Cd_rad_out_singleH;Cs_H_out_H/OneDe]
Euclidian distance = 1.41421356237
Multiplied by reaction path degeneracy 2.0
family: intra_H_migration"""),
)

reaction(
    label = 'reaction68',
    reactants = ['[CH]=[C]C([C]=C)C=O(24973)'],
    products = ['[CH]C(C=O)=C[C]=C(22615)'],
    transitionState = 'TS68',
    kinetics = Arrhenius(A=(8.889e+11,'s^-1'), n=0.232, Ea=(122.75,'kJ/mol'), T0=(1,'K'), comment="""Estimated using an average for rate rule [cCsCJ;CJ;CO]
Euclidian distance = 0
family: 1,2_shiftC"""),
)

reaction(
    label = 'reaction69',
    reactants = ['[CH]C(C=O)=C[C]=C(22615)'],
    products = ['C=[C]C1C=C1C=O(24974)'],
    transitionState = 'TS69',
    kinetics = Arrhenius(A=(2.18e+16,'s^-1'), n=0, Ea=(2.9288,'kJ/mol'), T0=(1,'K'), comment="""Estimated using template [R3;Y_rad_out;Ypri_rad_out] for rate rule [R3_SD;Y_rad_out;CdsinglepriH_rad_out]
Euclidian distance = 2.2360679775
family: Birad_recombination"""),
)

network(
    label = '4918',
    isomers = [
        '[CH]C(C=O)=C[C]=C(22615)',
    ],
    reactants = [
        ('C#CC=O(21959)', 'C#C[CH2](17441)'),
    ],
    bathGas = {
        'N2': 0.25,
        'Ne': 0.25,
        'He': 0.25,
        'Ar': 0.25,
    },
)

pressureDependence(
    label = '4918',
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

