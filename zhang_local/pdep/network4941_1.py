species(
    label = '[CH]=C=C[CH][C]=C[O](22638)',
    structure = SMILES('[CH]=C=C[CH][C]=C[O]'),
    E0 = (656.042,'kJ/mol'),
    modes = [
        HarmonicOscillator(frequencies=([3120,650,792.5,1650,3025,407.5,1350,352.5,540,610,2055,2995,3025,975,1000,1300,1375,400,500,1630,1680,1685,370,180,180],'cm^-1')),
        HinderedRotor(inertia=(1.84007,'amu*angstrom^2'), symmetry=1, barrier=(42.3067,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(1.84024,'amu*angstrom^2'), symmetry=1, barrier=(42.3107,'kJ/mol'), semiclassical=False),
    ],
    spinMultiplicity = 5,
    opticalIsomers = 1,
    molecularWeight = (92.0954,'amu'),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[1.30798,0.0451968,-1.54147e-06,-4.51244e-08,2.53584e-11,79013.7,25.3872], Tmin=(100,'K'), Tmax=(902.276,'K')), NASAPolynomial(coeffs=[18.8022,0.00563787,1.05541e-06,-3.69868e-10,2.58024e-14,74310.1,-65.784], Tmin=(902.276,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(656.042,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(249.434,'J/(mol*K)'), comment="""Thermo group additivity estimation: group(O2s-(Cds-Cd)H) + group(Cs-(Cds-Cds)(Cds-Cds)HH) + group(Cds-CdsCsH) + group(Cds-CdsCsH) + group(Cds-CdsOsH) + group(Cds-CdsHH) + group(Cdd-CdsCds) + radical(C=COJ) + radical(C=C=CJ) + radical(C=CCJC=C) + radical(Cds_S)"""),
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
    label = '[CH]=[C]C1[CH]C1=C[O](26150)',
    structure = SMILES('[CH]=[C]C1[CH]C1=C[O]'),
    E0 = (817.367,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
    molecularWeight = (92.0954,'amu'),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[1.18567,0.0510766,-2.5453e-05,-1.28943e-08,1.11889e-11,98417.7,20.9558], Tmin=(100,'K'), Tmax=(949.223,'K')), NASAPolynomial(coeffs=[17.0718,0.0104704,-2.90513e-06,5.00377e-10,-3.75041e-14,94215.3,-61.1078], Tmin=(949.223,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(817.367,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(253.591,'J/(mol*K)'), comment="""Thermo group additivity estimation: group(O2s-(Cds-Cd)H) + group(Cs-(Cds-Cds)(Cds-Cds)CsH) + group(Cs-(Cds-Cds)CsHH) + group(Cds-CdsCsCs) + group(Cds-CdsCsH) + group(Cds-CdsOsH) + group(Cds-CdsHH) + ring(Methylene_cyclopropane) + radical(Cds_S) + radical(Allyl_S) + radical(Cds_P) + radical(C=COJ)"""),
)

species(
    label = '[CH]=[C]C1[CH][C]=CO1(25321)',
    structure = SMILES('[CH]=[C]C1[CH][C]=CO1'),
    E0 = (731.642,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
    molecularWeight = (92.0954,'amu'),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[1.87476,0.0187553,9.00275e-05,-1.51378e-07,6.59188e-11,88098.9,20.9517], Tmin=(100,'K'), Tmax=(907.589,'K')), NASAPolynomial(coeffs=[23.8401,-0.00464489,7.37925e-06,-1.55179e-09,1.00609e-13,81088.5,-99.5444], Tmin=(907.589,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(731.642,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(253.591,'J/(mol*K)'), comment="""Thermo group additivity estimation: group(O2s-Cs(Cds-Cd)) + group(Cs-(Cds-Cds)CsOsH) + group(Cs-(Cds-Cds)CsHH) + group(Cds-CdsCsH) + group(Cds-CdsCsH) + group(Cds-CdsOsH) + group(Cds-CdsHH) + ring(2,3-Dihydrofuran) + radical(Cds_S) + radical(C=CCJC(O)C=C) + radical(Cds_P) + radical(Cds_S)"""),
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
    label = '[CH]=C=C=C[C]=C[O](26164)',
    structure = SMILES('[CH]=C=C=C[C]=C[O]'),
    E0 = (661.221,'kJ/mol'),
    modes = [
        HarmonicOscillator(frequencies=([3120,650,792.5,1650,1685,370,2995,3025,975,1000,1300,1375,400,500,1630,1680,540,563.333,586.667,610,1970,2140,180],'cm^-1')),
        HinderedRotor(inertia=(1.71577,'amu*angstrom^2'), symmetry=1, barrier=(39.4489,'kJ/mol'), semiclassical=False),
    ],
    spinMultiplicity = 4,
    opticalIsomers = 1,
    molecularWeight = (91.0874,'amu'),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[0.206843,0.0705449,-8.45013e-05,4.7546e-08,-9.83939e-12,79674.2,24.813], Tmin=(100,'K'), Tmax=(1380.72,'K')), NASAPolynomial(coeffs=[19.0489,0.00358679,1.68208e-06,-5.56426e-10,4.53257e-14,75650.4,-67.8967], Tmin=(1380.72,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(661.221,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(228.648,'J/(mol*K)'), comment="""Thermo group additivity estimation: group(O2s-(Cds-Cd)H) + group(Cds-Cds(Cds-Cds)H) + group(Cds-Cds(Cds-Cds)H) + group(Cds-CdsOsH) + group(Cds-CdsHH) + group(Cdd-CdsCds) + group(Cdd-CdsCds) + radical(C=C=CJ) + radical(C=CJC=C) + radical(C=COJ)"""),
)

species(
    label = '[CH]=C=C[CH][C]=C=O(26165)',
    structure = SMILES('[CH]=C=C[CH][C]=C=O'),
    E0 = (600.733,'kJ/mol'),
    modes = [
        HarmonicOscillator(frequencies=([3120,650,792.5,1650,3025,407.5,1350,352.5,540,610,2055,2120,512.5,787.5,1685,370,3010,987.5,1337.5,450,1655,180],'cm^-1')),
        HinderedRotor(inertia=(1.86658,'amu*angstrom^2'), symmetry=1, barrier=(42.9164,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(1.87259,'amu*angstrom^2'), symmetry=1, barrier=(43.0545,'kJ/mol'), semiclassical=False),
    ],
    spinMultiplicity = 4,
    opticalIsomers = 1,
    molecularWeight = (91.0874,'amu'),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[1.12351,0.0592714,-6.6441e-05,3.79538e-08,-8.44493e-12,72358.5,22.1235], Tmin=(100,'K'), Tmax=(1107.48,'K')), NASAPolynomial(coeffs=[13.7543,0.0136517,-4.65269e-06,7.5942e-10,-4.87867e-14,69560.8,-40.1022], Tmin=(1107.48,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(600.733,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(224.491,'J/(mol*K)'), comment="""Thermo group additivity estimation: group(Cs-Cd(CCO)HH) + group(Cds-CdsCsH) + group(Cds-(Cdd-O2d)CsH) + group(Cds-CdsHH) + group(Cdd-CdsCds) + radical(C=CCJC=C=O) + radical(C=C=CJ) + radical(CCCJ=C=O)"""),
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
    label = '[CH]=C=CC=C=[CH](21250)',
    structure = SMILES('[CH]=C=CC=C=[CH]'),
    E0 = (684.032,'kJ/mol'),
    modes = [
        HarmonicOscillator(frequencies=([3115,3125,620,680,785,800,1600,1700,3010,987.5,1337.5,450,1655,540,563.333,586.667,610,1970,2140,180,180,180,632.316],'cm^-1')),
        HinderedRotor(inertia=(4.39709,'amu*angstrom^2'), symmetry=1, barrier=(101.098,'kJ/mol'), semiclassical=False),
    ],
    spinMultiplicity = 3,
    opticalIsomers = 1,
    molecularWeight = (76.096,'amu'),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[0.96578,0.0556582,-6.00573e-05,3.22164e-08,-6.44519e-12,82388.8,21.1595], Tmin=(100,'K'), Tmax=(1428.06,'K')), NASAPolynomial(coeffs=[14.5868,0.0082227,-4.81583e-07,-1.47515e-10,1.73989e-14,79445.1,-46.093], Tmin=(1428.06,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(684.032,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(228.648,'J/(mol*K)'), comment="""Thermo group additivity estimation: group(Cds-Cds(Cds-Cds)H) + group(Cds-Cds(Cds-Cds)H) + group(Cds-CdsHH) + group(Cds-CdsHH) + group(Cdd-CdsCds) + group(Cdd-CdsCds) + radical(C=C=CJ) + radical(C=C=CJ)"""),
)

species(
    label = '[CH]=C=[C]C[C]=C[O](26166)',
    structure = SMILES('[CH]=C=[C]C[C]=C[O]'),
    E0 = (792.159,'kJ/mol'),
    modes = [
        HarmonicOscillator(frequencies=([3120,650,792.5,1650,2750,2850,1437.5,1250,1305,750,350,540,610,2055,3010,987.5,1337.5,450,1655,1670,1700,300,440,180,180],'cm^-1')),
        HinderedRotor(inertia=(1.12561,'amu*angstrom^2'), symmetry=1, barrier=(25.88,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(1.12437,'amu*angstrom^2'), symmetry=1, barrier=(25.8514,'kJ/mol'), semiclassical=False),
    ],
    spinMultiplicity = 5,
    opticalIsomers = 1,
    molecularWeight = (92.0954,'amu'),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[0.744706,0.0647338,-7.31319e-05,4.11705e-08,-8.88812e-12,95397.8,28.2082], Tmin=(100,'K'), Tmax=(1199.27,'K')), NASAPolynomial(coeffs=[16.0554,0.0111418,-2.94265e-06,3.96821e-10,-2.24176e-14,91907.1,-47.6819], Tmin=(1199.27,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(792.159,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(249.434,'J/(mol*K)'), comment="""Thermo group additivity estimation: group(O2s-(Cds-Cd)H) + group(Cs-(Cds-Cds)(Cds-Cds)HH) + group(Cds-CdsCsH) + group(Cds-CdsCsH) + group(Cds-CdsOsH) + group(Cds-CdsHH) + group(Cdd-CdsCds) + radical(Cds_S) + radical(Cds_S) + radical(C=C=CJ) + radical(C=COJ)"""),
)

species(
    label = '[CH]=C=C[CH][C]=[C]O(26167)',
    structure = SMILES('[CH]=C=C[CH][C]=[C]O'),
    E0 = (754.323,'kJ/mol'),
    modes = [
        HarmonicOscillator(frequencies=([3120,650,792.5,1650,3615,1277.5,1000,3025,407.5,1350,352.5,540,610,2055,1670,1700,300,440,3010,987.5,1337.5,450,1655,180],'cm^-1')),
        HinderedRotor(inertia=(1.23194,'amu*angstrom^2'), symmetry=1, barrier=(28.3248,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(1.80593,'amu*angstrom^2'), symmetry=1, barrier=(41.5218,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(1.22614,'amu*angstrom^2'), symmetry=1, barrier=(28.1914,'kJ/mol'), semiclassical=False),
    ],
    spinMultiplicity = 5,
    opticalIsomers = 1,
    molecularWeight = (92.0954,'amu'),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[0.231888,0.0642206,-6.7888e-05,3.45354e-08,-6.45164e-12,90875.9,30.9691], Tmin=(100,'K'), Tmax=(1565.19,'K')), NASAPolynomial(coeffs=[17.4451,0.0063574,8.60438e-07,-4.09787e-10,3.44922e-14,87186.8,-54.3578], Tmin=(1565.19,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(754.323,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(245.277,'J/(mol*K)'), comment="""Thermo group additivity estimation: group(O2s-(Cds-Cd)H) + group(Cs-(Cds-Cds)(Cds-Cds)HH) + group(Cds-CdsCsH) + group(Cds-CdsCsH) + group(Cds-CdsOsH) + group(Cds-CdsHH) + group(Cdd-CdsCds) + radical(C=CCJC=C) + radical(Cds_S) + radical(C=C=CJ) + radical(C=CJO)"""),
)

species(
    label = '[CH]=C=C[CH]C=[C][O](26168)',
    structure = SMILES('[CH]=C=C[CH]C=[C][O]'),
    E0 = (657.944,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
    molecularWeight = (92.0954,'amu'),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[1.5489,0.0423857,-3.78835e-06,-3.53293e-08,1.99365e-11,79231.4,27.1769], Tmin=(100,'K'), Tmax=(909.849,'K')), NASAPolynomial(coeffs=[15.8425,0.0102788,-1.52256e-06,1.34906e-10,-8.76521e-15,75358.4,-47.4216], Tmin=(909.849,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(657.944,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(249.434,'J/(mol*K)'), comment="""Thermo group additivity estimation: group(O2s-(Cds-Cd)H) + group(Cs-(Cds-Cds)(Cds-Cds)HH) + group(Cds-CdsCsH) + group(Cds-CdsCsH) + group(Cds-CdsOsH) + group(Cds-CdsHH) + group(Cdd-CdsCds) + radical(C=C=CJ) + radical(C=CCJC=C) + radical(C=COJ) + radical(C=CJO)"""),
)

species(
    label = '[CH]=C=[C][CH]C=C[O](26169)',
    structure = SMILES('[CH]=C=[C][CH]C=C[O]'),
    E0 = (656.042,'kJ/mol'),
    modes = [
        HarmonicOscillator(frequencies=([3120,650,792.5,1650,3025,407.5,1350,352.5,540,610,2055,2995,3025,975,1000,1300,1375,400,500,1630,1680,1685,370,180,180],'cm^-1')),
        HinderedRotor(inertia=(1.84007,'amu*angstrom^2'), symmetry=1, barrier=(42.3067,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(1.84024,'amu*angstrom^2'), symmetry=1, barrier=(42.3107,'kJ/mol'), semiclassical=False),
    ],
    spinMultiplicity = 5,
    opticalIsomers = 1,
    molecularWeight = (92.0954,'amu'),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[1.30798,0.0451968,-1.54147e-06,-4.51244e-08,2.53584e-11,79013.7,25.3872], Tmin=(100,'K'), Tmax=(902.276,'K')), NASAPolynomial(coeffs=[18.8022,0.00563787,1.05541e-06,-3.69868e-10,2.58024e-14,74310.1,-65.784], Tmin=(902.276,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(656.042,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(249.434,'J/(mol*K)'), comment="""Thermo group additivity estimation: group(O2s-(Cds-Cd)H) + group(Cs-(Cds-Cds)(Cds-Cds)HH) + group(Cds-CdsCsH) + group(Cds-CdsCsH) + group(Cds-CdsOsH) + group(Cds-CdsHH) + group(Cdd-CdsCds) + radical(C=C=CJ) + radical(C=COJ) + radical(Cds_S) + radical(C=CCJC=C)"""),
)

species(
    label = 'C#C[CH]C[C][C]=O(26170)',
    structure = SMILES('C#C[CH]C[C][C]=O'),
    E0 = (773.608,'kJ/mol'),
    modes = [
        HarmonicOscillator(frequencies=([1855,455,950,2750,2850,1437.5,1250,1305,750,350,750,770,3400,2100,2175,525,3025,407.5,1350,352.5,218.967,219.551,2750.57],'cm^-1')),
        HinderedRotor(inertia=(2.1308,'amu*angstrom^2'), symmetry=1, barrier=(72.3251,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(2.11164,'amu*angstrom^2'), symmetry=1, barrier=(72.3325,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.358669,'amu*angstrom^2'), symmetry=1, barrier=(12.2778,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(2.12058,'amu*angstrom^2'), symmetry=1, barrier=(72.332,'kJ/mol'), semiclassical=False),
    ],
    spinMultiplicity = 5,
    opticalIsomers = 1,
    molecularWeight = (92.0954,'amu'),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[0.901249,0.0723792,-0.000109962,8.7239e-08,-2.64387e-11,93151.2,26.126], Tmin=(100,'K'), Tmax=(947.422,'K')), NASAPolynomial(coeffs=[10.3573,0.020013,-7.35372e-06,1.17484e-09,-7.07252e-14,91917.8,-16.036], Tmin=(947.422,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(773.608,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(241.12,'J/(mol*K)'), comment="""Thermo group additivity estimation: group(Cs-CsCsHH) + group(Cs-(Cds-O2d)CsHH) + group(Cs-CtCsHH) + group(Cds-OdCsH) + group(Ct-CtCs) + group(Ct-CtH) + radical(Sec_Propargyl) + radical(CCCJ=O) + radical(CCJ2_triplet)"""),
)

species(
    label = 'C=[C][C]=C[C]=C[O](24988)',
    structure = SMILES('C=[C][C]=C[C]=C[O]'),
    E0 = (675.35,'kJ/mol'),
    modes = [
        HarmonicOscillator(frequencies=([2950,3100,1380,975,1025,1650,2995,3025,975,1000,1300,1375,400,500,1630,1680,1670,1685,1700,300,370,440,180,180,180],'cm^-1')),
        HinderedRotor(inertia=(1.77597,'amu*angstrom^2'), symmetry=1, barrier=(40.8331,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(1.78786,'amu*angstrom^2'), symmetry=1, barrier=(41.1063,'kJ/mol'), semiclassical=False),
    ],
    spinMultiplicity = 5,
    opticalIsomers = 1,
    molecularWeight = (92.0954,'amu'),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[-0.000283518,0.0762297,-9.20979e-05,5.33849e-08,-1.14246e-11,81380,25.0358], Tmin=(100,'K'), Tmax=(1336.68,'K')), NASAPolynomial(coeffs=[18.8718,0.00648097,1.06916e-06,-5.11427e-10,4.51732e-14,77520.7,-67.0517], Tmin=(1336.68,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(675.35,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(249.434,'J/(mol*K)'), comment="""Thermo group additivity estimation: group(O2s-(Cds-Cd)H) + group(Cds-Cds(Cds-Cds)H) + group(Cds-Cds(Cds-Cds)H) + group(Cds-Cds(Cds-Cds)H) + group(Cds-Cds(Cds-Cds)H) + group(Cds-CdsOsH) + group(Cds-CdsHH) + radical(C=CJC=C) + radical(C=CJC=C) + radical(C=CJC=C) + radical(C=COJ)"""),
)

species(
    label = '[CH]=C=[C][CH][C]=CO(26171)',
    structure = SMILES('[CH]=C=[C][CH][C]=CO'),
    E0 = (752.421,'kJ/mol'),
    modes = [
        HarmonicOscillator(frequencies=([3120,650,792.5,1650,3615,1277.5,1000,3025,407.5,1350,352.5,540,610,2055,1670,1700,300,440,3010,987.5,1337.5,450,1655,180],'cm^-1')),
        HinderedRotor(inertia=(1.45235,'amu*angstrom^2'), symmetry=1, barrier=(33.3923,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(1.45416,'amu*angstrom^2'), symmetry=1, barrier=(33.434,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(1.44796,'amu*angstrom^2'), symmetry=1, barrier=(33.2915,'kJ/mol'), semiclassical=False),
    ],
    spinMultiplicity = 5,
    opticalIsomers = 1,
    molecularWeight = (92.0954,'amu'),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[-0.287568,0.0702746,-7.68594e-05,3.92439e-08,-7.24817e-12,90670.4,30.1823], Tmin=(100,'K'), Tmax=(1611.32,'K')), NASAPolynomial(coeffs=[19.3882,0.00306955,2.79546e-06,-7.84376e-10,5.95249e-14,86713.3,-66.7319], Tmin=(1611.32,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(752.421,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(245.277,'J/(mol*K)'), comment="""Thermo group additivity estimation: group(O2s-(Cds-Cd)H) + group(Cs-(Cds-Cds)(Cds-Cds)HH) + group(Cds-CdsCsH) + group(Cds-CdsCsH) + group(Cds-CdsOsH) + group(Cds-CdsHH) + group(Cdd-CdsCds) + radical(C=CCJC=C) + radical(Cds_S) + radical(Cds_S) + radical(C=C=CJ)"""),
)

species(
    label = '[CH2][C]=C[CH][C]=C=O(24990)',
    structure = SMILES('[CH2][C]=C[CH][C]=C=O'),
    E0 = (658.993,'kJ/mol'),
    modes = [
        HarmonicOscillator(frequencies=([1670,1700,300,440,3025,407.5,1350,352.5,2120,512.5,787.5,3010,987.5,1337.5,450,1655,3000,3100,440,815,1455,1000,323.127,324.139],'cm^-1')),
        HinderedRotor(inertia=(0.546344,'amu*angstrom^2'), symmetry=1, barrier=(40.7222,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.545682,'amu*angstrom^2'), symmetry=1, barrier=(40.7354,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.54873,'amu*angstrom^2'), symmetry=1, barrier=(40.736,'kJ/mol'), semiclassical=False),
    ],
    spinMultiplicity = 5,
    opticalIsomers = 1,
    molecularWeight = (92.0954,'amu'),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[1.20327,0.0587517,-5.85757e-05,3.00532e-08,-6.14597e-12,79361.6,22.9137], Tmin=(100,'K'), Tmax=(1183.13,'K')), NASAPolynomial(coeffs=[13.0506,0.0186969,-7.79266e-06,1.43768e-09,-9.93312e-14,76558.2,-36.235], Tmin=(1183.13,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(658.993,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(245.277,'J/(mol*K)'), comment="""Thermo group additivity estimation: group(Cs-Cd(CCO)HH) + group(Cs-(Cds-Cds)HHH) + group(Cds-CdsCsH) + group(Cds-CdsCsH) + group(Cds-(Cdd-O2d)CsH) + radical(CCCJ=C=O) + radical(C=CCJC=C=O) + radical(Cds_S) + radical(Allyl_P)"""),
)

species(
    label = '[CH]=C=[C][CH][C]=C[O](26172)',
    structure = SMILES('[CH]=C=[C][CH][C]=C[O]'),
    E0 = (893.884,'kJ/mol'),
    modes = [
        HarmonicOscillator(frequencies=([3120,650,792.5,1650,3025,407.5,1350,352.5,540,610,2055,3010,987.5,1337.5,450,1655,1670,1700,300,440,180,180],'cm^-1')),
        HinderedRotor(inertia=(1.86111,'amu*angstrom^2'), symmetry=1, barrier=(42.7906,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(1.86897,'amu*angstrom^2'), symmetry=1, barrier=(42.9713,'kJ/mol'), semiclassical=False),
    ],
    spinMultiplicity = 6,
    opticalIsomers = 1,
    molecularWeight = (91.0874,'amu'),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[0.378206,0.0602626,-6.25834e-05,3.09671e-08,-5.62006e-12,107656,29.2057], Tmin=(100,'K'), Tmax=(1613.59,'K')), NASAPolynomial(coeffs=[17.0885,0.00545975,7.98779e-07,-3.5824e-10,2.94509e-14,104005,-54.0101], Tmin=(1613.59,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(893.884,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(224.491,'J/(mol*K)'), comment="""Thermo group additivity estimation: group(O2s-(Cds-Cd)H) + group(Cs-(Cds-Cds)(Cds-Cds)HH) + group(Cds-CdsCsH) + group(Cds-CdsCsH) + group(Cds-CdsOsH) + group(Cds-CdsHH) + group(Cdd-CdsCds) + radical(Cds_S) + radical(C=C=CJ) + radical(C=CCJC=C) + radical(C=COJ) + radical(Cds_S)"""),
)

species(
    label = '[CH][C]=C[CH][C]=C=O(26173)',
    structure = SMILES('[CH][C]=C[CH][C]=C=O'),
    E0 = (878.178,'kJ/mol'),
    modes = [
        HarmonicOscillator(frequencies=([3025,407.5,1350,352.5,2120,512.5,787.5,3010,987.5,1337.5,450,1655,1670,1700,300,440,256.391,257.007,257.598,258.621,260.871],'cm^-1')),
        HinderedRotor(inertia=(1.06277,'amu*angstrom^2'), symmetry=1, barrier=(50.3132,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(1.06085,'amu*angstrom^2'), symmetry=1, barrier=(50.3049,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(1.08244,'amu*angstrom^2'), symmetry=1, barrier=(50.3644,'kJ/mol'), semiclassical=False),
    ],
    spinMultiplicity = 6,
    opticalIsomers = 1,
    molecularWeight = (91.0874,'amu'),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[1.37243,0.0591149,-5.93078e-05,3.29345e-08,-7.59219e-12,105714,23.1399], Tmin=(100,'K'), Tmax=(1033.27,'K')), NASAPolynomial(coeffs=[9.85409,0.0262804,-1.16417e-05,2.18012e-09,-1.511e-13,103961,-18.0567], Tmin=(1033.27,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(878.178,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(220.334,'J/(mol*K)'), comment="""Thermo group additivity estimation: group(Cs-Cd(CCO)HH) + group(Cs-(Cds-Cds)HHH) + group(Cds-CdsCsH) + group(Cds-CdsCsH) + group(Cds-(Cdd-O2d)CsH) + radical(AllylJ2_triplet) + radical(CCCJ=C=O) + radical(Cds_S) + radical(C=CCJC=C=O)"""),
)

species(
    label = '[CH]=C1[CH]C1[C]=C[O](26085)',
    structure = SMILES('[CH]=C1[CH]C1[C]=C[O]'),
    E0 = (817.367,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
    molecularWeight = (92.0954,'amu'),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[1.18571,0.0510761,-2.54513e-05,-1.28966e-08,1.11899e-11,98417.7,20.9557], Tmin=(100,'K'), Tmax=(949.214,'K')), NASAPolynomial(coeffs=[17.0717,0.0104706,-2.90525e-06,5.00407e-10,-3.75066e-14,94215.3,-61.1071], Tmin=(949.214,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(817.367,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(253.591,'J/(mol*K)'), comment="""Thermo group additivity estimation: group(O2s-(Cds-Cd)H) + group(Cs-(Cds-Cds)(Cds-Cds)CsH) + group(Cs-(Cds-Cds)CsHH) + group(Cds-CdsCsCs) + group(Cds-CdsCsH) + group(Cds-CdsOsH) + group(Cds-CdsHH) + ring(Methylene_cyclopropane) + radical(Cds_S) + radical(C=COJ) + radical(Allyl_S) + radical(Cds_P)"""),
)

species(
    label = '[O]C=[C][CH]C1[C]=C1(26174)',
    structure = SMILES('[O]C=[C][CH]C1[C]=C1'),
    E0 = (861.146,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
    molecularWeight = (92.0954,'amu'),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[1.05995,0.0577335,-5.18627e-05,1.95482e-08,-1.62086e-12,103684,22.2732], Tmin=(100,'K'), Tmax=(1007.6,'K')), NASAPolynomial(coeffs=[15.3319,0.0136457,-4.94152e-06,8.83757e-10,-6.17185e-14,100170,-49.8549], Tmin=(1007.6,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(861.146,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(249.434,'J/(mol*K)'), comment="""Thermo group additivity estimation: group(O2s-(Cds-Cd)H) + group(Cs-(Cds-Cds)(Cds-Cds)CsH) + group(Cs-(Cds-Cds)CsHH) + group(Cds-CdsCsH) + group(Cds-CdsCsH) + group(Cds-CdsCsH) + group(Cds-CdsOsH) + ring(Cyclopropene) + radical(Allyl_S) + radical(C=COJ) + radical(Cds_S) + radical(cyclopropenyl-vinyl)"""),
)

species(
    label = '[CH]C1=C[CH]C1=C[O](26175)',
    structure = SMILES('[CH]C1=C[CH]C1=C[O]'),
    E0 = (587.213,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
    molecularWeight = (92.0954,'amu'),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[1.51848,0.025457,8.91497e-05,-1.53132e-07,6.63255e-11,70741.9,20.0706], Tmin=(100,'K'), Tmax=(915.113,'K')), NASAPolynomial(coeffs=[24.3207,0.00148826,4.35381e-06,-9.61136e-10,5.86201e-14,63398.9,-105.233], Tmin=(915.113,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(587.213,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(253.591,'J/(mol*K)'), comment="""Thermo group additivity estimation: group(O2s-(Cds-Cd)H) + group(Cs-(Cds-Cds)(Cds-Cds)HH) + group(Cs-(Cds-Cds)HHH) + group(Cds-Cds(Cds-Cds)Cs) + group(Cds-Cds(Cds-Cds)Cs) + group(Cds-CdsCsH) + group(Cds-CdsOsH) + ring(Cyclobutene) + radical(C=COJ) + radical(AllylJ2_triplet) + radical(C=CCJC=C)"""),
)

species(
    label = '[CH]C1=CC=[C][CH]O1(26110)',
    structure = SMILES('[CH]C1=CC=[C][CH]O1'),
    E0 = (626.398,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
    molecularWeight = (92.0954,'amu'),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[1.75471,0.0342876,2.81877e-05,-6.19271e-08,2.61253e-11,75432.7,19.29], Tmin=(100,'K'), Tmax=(973.822,'K')), NASAPolynomial(coeffs=[13.9905,0.0208199,-7.73802e-06,1.46299e-09,-1.07901e-13,71305.1,-48.3731], Tmin=(973.822,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(626.398,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(253.591,'J/(mol*K)'), comment="""Thermo group additivity estimation: group(O2s-Cs(Cds-Cd)) + group(Cs-(Cds-Cds)OsHH) + group(Cs-(Cds-Cds)HHH) + group(Cds-CdsCsOs) + group(Cds-CdsCsH) + group(Cds-Cds(Cds-Cds)H) + group(Cds-Cds(Cds-Cds)H) + ring(1,3-Cyclohexadiene) + radical(Cds_S) + radical(C=CCJ(O)C) + radical(AllylJ2_triplet)"""),
)

species(
    label = '[CH]=C=C[CH]C=C=O(26176)',
    structure = SMILES('[CH]=C=C[CH]C=C=O'),
    E0 = (398.452,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
    molecularWeight = (92.0954,'amu'),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[1.11708,0.0584635,-5.99207e-05,3.18861e-08,-6.67616e-12,48030.7,22.4381], Tmin=(100,'K'), Tmax=(1169.62,'K')), NASAPolynomial(coeffs=[13.4973,0.0161244,-5.62267e-06,9.37256e-10,-6.10616e-14,45134.7,-39.2293], Tmin=(1169.62,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(398.452,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(249.434,'J/(mol*K)'), comment="""Thermo group additivity estimation: group(Cs-Cd(CCO)HH) + group(Cds-CdsCsH) + group(Cds-(Cdd-O2d)CsH) + group(Cds-CdsHH) + group(Cdd-CdsCds) + radical(C=C=CJ) + radical(C=CCJC=C=O)"""),
)

species(
    label = '[CH2]C#CC=C=C[O](24976)',
    structure = SMILES('[CH2]C#CC=C=C[O]'),
    E0 = (458.182,'kJ/mol'),
    modes = [
        HarmonicOscillator(frequencies=([540,610,2055,2995,3025,975,1000,1300,1375,400,500,1630,1680,2100,2250,500,550,3000,3100,440,815,1455,1000,180,180],'cm^-1')),
        HinderedRotor(inertia=(1.57299,'amu*angstrom^2'), symmetry=1, barrier=(36.166,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(1.58402,'amu*angstrom^2'), symmetry=1, barrier=(36.4198,'kJ/mol'), semiclassical=False),
    ],
    spinMultiplicity = 3,
    opticalIsomers = 1,
    molecularWeight = (92.0954,'amu'),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[1.20454,0.0485282,-1.46144e-05,-2.56157e-08,1.581e-11,55219,23.5289], Tmin=(100,'K'), Tmax=(955.365,'K')), NASAPolynomial(coeffs=[18.0513,0.00944701,-2.63911e-06,4.89588e-10,-3.92024e-14,50564.6,-64.4902], Tmin=(955.365,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(458.182,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(249.434,'J/(mol*K)'), comment="""Thermo group additivity estimation: group(O2s-(Cds-Cd)H) + group(Cs-CtHHH) + group(Cds-CdsCtH) + group(Cds-CdsOsH) + group(Cdd-CdsCds) + group(Ct-CtCs) + group(Ct-Ct(Cds-Cds)) + radical(Propargyl) + radical(C=COJ)"""),
)

species(
    label = '[CH]=C=CC=C1[CH]O1(26177)',
    structure = SMILES('[CH]=C=CC=C1[CH]O1'),
    E0 = (488.407,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
    molecularWeight = (92.0954,'amu'),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[1.21279,0.0366731,4.41067e-05,-1.08056e-07,5.14607e-11,58865.3,20.1977], Tmin=(100,'K'), Tmax=(906.079,'K')), NASAPolynomial(coeffs=[25.8644,-0.00674552,7.70203e-06,-1.59863e-09,1.05024e-13,51713.1,-111.117], Tmin=(906.079,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(488.407,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(253.591,'J/(mol*K)'), comment="""Thermo group additivity estimation: group(O2s-Cs(Cds-Cd)) + group(Cs-(Cds-Cds)OsHH) + group(Cds-CdsCsOs) + group(Cds-Cds(Cds-Cds)H) + group(Cds-Cds(Cds-Cds)H) + group(Cds-CdsHH) + group(Cdd-CdsCds) + ring(methyleneoxirane) + radical(C=CCJO) + radical(C=C=CJ)"""),
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
    label = '[CH][C]=CC=C=[CH](21252)',
    structure = SMILES('[CH][C]=CC=C=[CH]'),
    E0 = (961.477,'kJ/mol'),
    modes = [
        HarmonicOscillator(frequencies=([3120,650,792.5,1650,1685,370,2995,3025,975,1000,1300,1375,400,500,1630,1680,540,610,2055,180,180,180],'cm^-1')),
        HinderedRotor(inertia=(2.0902,'amu*angstrom^2'), symmetry=1, barrier=(48.0578,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(2.08995,'amu*angstrom^2'), symmetry=1, barrier=(48.052,'kJ/mol'), semiclassical=False),
    ],
    spinMultiplicity = 5,
    opticalIsomers = 1,
    molecularWeight = (76.096,'amu'),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[1.38968,0.053466,-4.58862e-05,1.80981e-08,-1.6973e-12,115737,21.5458], Tmin=(100,'K'), Tmax=(934.723,'K')), NASAPolynomial(coeffs=[11.5298,0.019595,-6.81242e-06,1.1284e-09,-7.35358e-14,113425,-28.9152], Tmin=(934.723,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(961.477,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(224.491,'J/(mol*K)'), comment="""Thermo group additivity estimation: group(Cs-(Cds-Cds)HHH) + group(Cds-CdsCsH) + group(Cds-Cds(Cds-Cds)H) + group(Cds-Cds(Cds-Cds)H) + group(Cds-CdsHH) + group(Cdd-CdsCds) + radical(C=C=CJ) + radical(Cds_S) + radical(AllylJ2_triplet)"""),
)

species(
    label = '[CH]=[C]C=[C][C]=C[O](26178)',
    structure = SMILES('[CH]=[C]C=[C][C]=C[O]'),
    E0 = (922.446,'kJ/mol'),
    modes = [
        HarmonicOscillator(frequencies=([1670,1685,1700,300,370,440,2995,3025,975,1000,1300,1375,400,500,1630,1680,3120,650,792.5,1650,180,180],'cm^-1')),
        HinderedRotor(inertia=(1.76796,'amu*angstrom^2'), symmetry=1, barrier=(40.6488,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(1.76786,'amu*angstrom^2'), symmetry=1, barrier=(40.6466,'kJ/mol'), semiclassical=False),
    ],
    spinMultiplicity = 6,
    opticalIsomers = 1,
    molecularWeight = (91.0874,'amu'),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[0.0568389,0.0782343,-0.000103568,6.42215e-08,-1.46326e-11,111094,25.2325], Tmin=(100,'K'), Tmax=(1253.59,'K')), NASAPolynomial(coeffs=[19.0912,0.00374197,2.02759e-06,-6.88942e-10,5.80801e-14,107403,-66.588], Tmin=(1253.59,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(922.446,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(224.491,'J/(mol*K)'), comment="""Thermo group additivity estimation: group(O2s-(Cds-Cd)H) + group(Cds-Cds(Cds-Cds)H) + group(Cds-Cds(Cds-Cds)H) + group(Cds-Cds(Cds-Cds)H) + group(Cds-Cds(Cds-Cds)H) + group(Cds-CdsOsH) + group(Cds-CdsHH) + radical(Cds_P) + radical(C=CJC=C) + radical(C=COJ) + radical(C=CJC=C) + radical(C=CJC=C)"""),
)

species(
    label = '[C]#C[CH][CH][C]=C[O](26179)',
    structure = SMILES('[C]#C[CH][CH][C]=C[O]'),
    E0 = (1024.78,'kJ/mol'),
    modes = [
        HarmonicOscillator(frequencies=([3000,3050,390,425,1340,1360,335,370,3010,987.5,1337.5,450,1655,2175,525,1685,370,400.853,400.855,400.857,400.86],'cm^-1')),
        HinderedRotor(inertia=(0.489312,'amu*angstrom^2'), symmetry=1, barrier=(55.795,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.489326,'amu*angstrom^2'), symmetry=1, barrier=(55.7952,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.00104913,'amu*angstrom^2'), symmetry=1, barrier=(0.119627,'kJ/mol'), semiclassical=False),
    ],
    spinMultiplicity = 6,
    opticalIsomers = 1,
    molecularWeight = (91.0874,'amu'),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[0.844947,0.0625034,-7.55337e-05,4.53746e-08,-1.0188e-11,123372,25.4149], Tmin=(100,'K'), Tmax=(1261.09,'K')), NASAPolynomial(coeffs=[14.6868,0.00910203,-7.19696e-07,-1.46764e-10,1.99811e-14,120636,-41.581], Tmin=(1261.09,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(1024.78,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(220.334,'J/(mol*K)'), comment="""Thermo group additivity estimation: group(O2s-(Cds-Cd)H) + group(Cs-(Cds-Cds)CsHH) + group(Cs-CtCsHH) + group(Cds-CdsCsH) + group(Cds-CdsOsH) + group(Ct-CtCs) + group(Ct-CtH) + radical(C=COJ) + radical(Sec_Propargyl) + radical(Allyl_S) + radical(Acetyl) + radical(Cds_S)"""),
)

species(
    label = '[C]=[C]CC=C=C[O](24806)',
    structure = SMILES('[C]=[C]CC=C=C[O]'),
    E0 = (957.942,'kJ/mol'),
    modes = [
        HarmonicOscillator(frequencies=([1685,370,2750,2850,1437.5,1250,1305,750,350,2995,3025,975,1000,1300,1375,400,500,1630,1680,540,610,2055,180,180,180],'cm^-1')),
        HinderedRotor(inertia=(0.994823,'amu*angstrom^2'), symmetry=1, barrier=(22.8729,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.999535,'amu*angstrom^2'), symmetry=1, barrier=(22.9813,'kJ/mol'), semiclassical=False),
    ],
    spinMultiplicity = 5,
    opticalIsomers = 1,
    molecularWeight = (92.0954,'amu'),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[0.954221,0.0627962,-6.78491e-05,3.64361e-08,-7.63102e-12,115327,26.7374], Tmin=(100,'K'), Tmax=(1170.65,'K')), NASAPolynomial(coeffs=[15.1503,0.0142891,-5.69467e-06,1.03992e-09,-7.18934e-14,112003,-43.9872], Tmin=(1170.65,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(957.942,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(249.434,'J/(mol*K)'), comment="""Thermo group additivity estimation: group(O2s-(Cds-Cd)H) + group(Cs-(Cds-Cds)(Cds-Cds)HH) + group(Cds-CdsCsH) + group(Cds-CdsCsH) + group(Cds-CdsHH) + group(Cds-CdsOsH) + group(Cdd-CdsCds) + radical(Cds_S) + radical(C=COJ) + radical(CdCdJ2_triplet)"""),
)

species(
    label = '[CH]=[C]C[CH][C]=C=O(24794)',
    structure = SMILES('[CH]=[C]C[CH][C]=C=O'),
    E0 = (813.059,'kJ/mol'),
    modes = [
        HarmonicOscillator(frequencies=([1670,1700,300,440,2750,2850,1437.5,1250,1305,750,350,2120,512.5,787.5,3120,650,792.5,1650,3025,407.5,1350,352.5,292.55,293.789],'cm^-1')),
        HinderedRotor(inertia=(0.0688635,'amu*angstrom^2'), symmetry=1, barrier=(4.11457,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.355583,'amu*angstrom^2'), symmetry=1, barrier=(21.6174,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.914852,'amu*angstrom^2'), symmetry=1, barrier=(55.5995,'kJ/mol'), semiclassical=False),
    ],
    spinMultiplicity = 5,
    opticalIsomers = 1,
    molecularWeight = (92.0954,'amu'),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[1.18871,0.0666615,-9.54696e-05,7.79543e-08,-2.59466e-11,97885.1,27.1924], Tmin=(100,'K'), Tmax=(771.998,'K')), NASAPolynomial(coeffs=[8.45671,0.0265588,-1.28001e-05,2.4625e-09,-1.71421e-13,96835.8,-5.51889], Tmin=(771.998,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(813.059,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(245.277,'J/(mol*K)'), comment="""Thermo group additivity estimation: group(Cs-(Cds-Cds)CsHH) + group(Cs-(Cds-Cdd-O2d)CsHH) + group(Cds-CdsCsH) + group(Cds-(Cdd-O2d)CsH) + group(Cds-CdsHH) + radical(CCCJ=C=O) + radical(Cds_S) + radical(CCJC(C)=C=O) + radical(Cds_P)"""),
)

species(
    label = '[C]#C[CH]C[C]=C[O](26180)',
    structure = SMILES('[C]#C[CH]C[C]=C[O]'),
    E0 = (883.669,'kJ/mol'),
    modes = [
        HarmonicOscillator(frequencies=([1685,370,2750,2850,1437.5,1250,1305,750,350,3010,987.5,1337.5,450,1655,2175,525,3025,407.5,1350,352.5,264.829,264.83,264.836,264.86],'cm^-1')),
        HinderedRotor(inertia=(1.25226,'amu*angstrom^2'), symmetry=1, barrier=(62.3286,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.00240333,'amu*angstrom^2'), symmetry=1, barrier=(0.119627,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(1.25202,'amu*angstrom^2'), symmetry=1, barrier=(62.3278,'kJ/mol'), semiclassical=False),
    ],
    spinMultiplicity = 5,
    opticalIsomers = 1,
    molecularWeight = (92.0954,'amu'),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[1.06611,0.0623983,-6.37998e-05,2.35184e-08,1.64554e-12,106389,25.752], Tmin=(100,'K'), Tmax=(792.966,'K')), NASAPolynomial(coeffs=[14.7133,0.0115843,-1.77942e-06,4.59914e-11,6.97278e-15,103658,-40.4959], Tmin=(792.966,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(883.669,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(245.277,'J/(mol*K)'), comment="""Thermo group additivity estimation: group(O2s-(Cds-Cd)H) + group(Cs-(Cds-Cds)CsHH) + group(Cs-CtCsHH) + group(Cds-CdsCsH) + group(Cds-CdsOsH) + group(Ct-CtCs) + group(Ct-CtH) + radical(Acetyl) + radical(Sec_Propargyl) + radical(C=COJ) + radical(Cds_S)"""),
)

species(
    label = '[C]#C[CH][CH]C=C[O](26181)',
    structure = SMILES('[C]#C[CH][CH]C=C[O]'),
    E0 = (786.94,'kJ/mol'),
    modes = [
        HarmonicOscillator(frequencies=([3000,3050,390,425,1340,1360,335,370,2995,3025,975,1000,1300,1375,400,500,1630,1680,2175,525,504.329,504.33,504.33,504.33],'cm^-1')),
        HinderedRotor(inertia=(0.382939,'amu*angstrom^2'), symmetry=1, barrier=(69.117,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.382939,'amu*angstrom^2'), symmetry=1, barrier=(69.117,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.382939,'amu*angstrom^2'), symmetry=1, barrier=(69.117,'kJ/mol'), semiclassical=False),
    ],
    spinMultiplicity = 5,
    opticalIsomers = 1,
    molecularWeight = (92.0954,'amu'),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[0.516195,0.0620753,-6.49193e-05,3.39317e-08,-6.59474e-12,94784.9,26.1272], Tmin=(100,'K'), Tmax=(1489.35,'K')), NASAPolynomial(coeffs=[15.5927,0.00984689,-4.96643e-07,-1.96533e-10,2.21658e-14,91595.8,-48.2441], Tmin=(1489.35,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(786.94,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(245.277,'J/(mol*K)'), comment="""Thermo group additivity estimation: group(O2s-(Cds-Cd)H) + group(Cs-(Cds-Cds)CsHH) + group(Cs-CtCsHH) + group(Cds-CdsCsH) + group(Cds-CdsOsH) + group(Ct-CtCs) + group(Ct-CtH) + radical(Acetyl) + radical(Sec_Propargyl) + radical(C=COJ) + radical(Allyl_S)"""),
)

species(
    label = '[C]#C[CH][CH][C]=CO(26182)',
    structure = SMILES('[C]#C[CH][CH][C]=CO'),
    E0 = (883.319,'kJ/mol'),
    modes = [
        HarmonicOscillator(frequencies=([1685,370,3000,3050,390,425,1340,1360,335,370,3615,1277.5,1000,2175,525,3010,987.5,1337.5,450,1655,380.931,380.931,380.931],'cm^-1')),
        HinderedRotor(inertia=(0.755701,'amu*angstrom^2'), symmetry=1, barrier=(77.816,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.158136,'amu*angstrom^2'), symmetry=1, barrier=(16.2836,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.755702,'amu*angstrom^2'), symmetry=1, barrier=(77.816,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.7557,'amu*angstrom^2'), symmetry=1, barrier=(77.816,'kJ/mol'), semiclassical=False),
    ],
    spinMultiplicity = 5,
    opticalIsomers = 1,
    molecularWeight = (92.0954,'amu'),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[0.215278,0.0721252,-8.8601e-05,5.22902e-08,-1.13194e-11,106385,26.2595], Tmin=(100,'K'), Tmax=(1340.83,'K')), NASAPolynomial(coeffs=[17.5299,0.00596394,1.64404e-06,-6.49478e-10,5.5818e-14,103046,-57.4888], Tmin=(1340.83,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(883.319,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(241.12,'J/(mol*K)'), comment="""Thermo group additivity estimation: group(O2s-(Cds-Cd)H) + group(Cs-(Cds-Cds)CsHH) + group(Cs-CtCsHH) + group(Cds-CdsCsH) + group(Cds-CdsOsH) + group(Ct-CtCs) + group(Ct-CtH) + radical(Acetyl) + radical(Sec_Propargyl) + radical(Cds_S) + radical(Allyl_S)"""),
)

species(
    label = '[O]C=[C]C1[CH][C]=C1(26183)',
    structure = SMILES('[O]C=[C]C1[CH][C]=C1'),
    E0 = (803.71,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
    molecularWeight = (92.0954,'amu'),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[1.57393,0.0412766,-3.02707e-06,-3.17198e-08,1.66731e-11,96762.3,20.0637], Tmin=(100,'K'), Tmax=(959.925,'K')), NASAPolynomial(coeffs=[15.477,0.0126231,-4.0068e-06,7.36913e-10,-5.56127e-14,92744,-53.4688], Tmin=(959.925,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(803.71,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(253.591,'J/(mol*K)'), comment="""Thermo group additivity estimation: group(O2s-(Cds-Cd)H) + group(Cs-(Cds-Cds)(Cds-Cds)CsH) + group(Cs-(Cds-Cds)CsHH) + group(Cds-CdsCsH) + group(Cds-CdsCsH) + group(Cds-CdsCsH) + group(Cds-CdsOsH) + ring(Cyclobutene) + radical(cyclobutene-vinyl) + radical(cyclobutene-allyl) + radical(Cds_S) + radical(C=COJ)"""),
)

species(
    label = '[O]C=C1[CH][C]=C[CH]1(26160)',
    structure = SMILES('[O]C=C1[CH][C]=C[CH]1'),
    E0 = (508.997,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
    molecularWeight = (92.0954,'amu'),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[2.59207,0.00100813,0.000131651,-1.84143e-07,7.43734e-11,61296.7,20.1401], Tmin=(100,'K'), Tmax=(922.031,'K')), NASAPolynomial(coeffs=[19.9994,0.00263682,3.49705e-06,-7.36716e-10,3.98696e-14,54807.4,-80.21], Tmin=(922.031,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(508.997,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(257.749,'J/(mol*K)'), comment="""Thermo group additivity estimation: group(O2s-(Cds-Cd)H) + group(Cs-(Cds-Cds)(Cds-Cds)HH) + group(Cs-(Cds-Cds)(Cds-Cds)HH) + group(Cds-CdsCsCs) + group(Cds-CdsCsH) + group(Cds-CdsCsH) + group(Cds-CdsOsH) + ring(4-Methylenecyclopentene) + radical(C=CCJC=C) + radical(C=CCJC=C) + radical(C=COJ) + radical(cyclopentene-vinyl)"""),
)

species(
    label = '[C]1[CH]OC=[C][CH]C=1(25334)',
    structure = SMILES('[C]1[CH]OC=[C][CH]C=1'),
    E0 = (647.834,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
    molecularWeight = (92.0954,'amu'),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[2.35512,-0.02763,0.000284009,-4.0218e-07,1.68867e-10,78035.9,20.5441], Tmin=(100,'K'), Tmax=(896.724,'K')), NASAPolynomial(coeffs=[43.6872,-0.0444477,3.18684e-05,-6.35901e-09,4.24214e-13,63886.7,-211.915], Tmin=(896.724,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(647.834,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(257.749,'J/(mol*K)'), comment="""Thermo group additivity estimation: group(O2s-Cs(Cds-Cd)) + group(Cs-(Cds-Cds)(Cds-Cds)HH) + group(Cs-(Cds-Cds)OsHH) + group(Cds-CdsCsH) + group(Cds-CdsCsH) + group(Cds-CdsCsH) + group(Cds-CdsOsH) + ring(Cycloheptane) + radical(C=CCJ(O)C) + radical(Cds_S) + radical(C=CCJC=C) + radical(Cds_S)"""),
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
    label = '[CH][CH][C]=C[O](26184)',
    structure = SMILES('[CH][CH][C]=C[O]'),
    E0 = (743.318,'kJ/mol'),
    modes = [
        HarmonicOscillator(frequencies=([3025,407.5,1350,352.5,3010,987.5,1337.5,450,1655,1685,370,446.686,446.948,447.724,447.75,448.081],'cm^-1')),
        HinderedRotor(inertia=(0.000839229,'amu*angstrom^2'), symmetry=1, barrier=(0.119627,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.00084471,'amu*angstrom^2'), symmetry=1, barrier=(0.119627,'kJ/mol'), semiclassical=False),
    ],
    spinMultiplicity = 6,
    opticalIsomers = 1,
    molecularWeight = (67.066,'amu'),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[2.14993,0.031358,-5.14644e-06,-2.37988e-08,1.3616e-11,89475.7,18.359], Tmin=(100,'K'), Tmax=(934.877,'K')), NASAPolynomial(coeffs=[13.8477,0.00495333,-7.20062e-07,9.96693e-11,-9.65933e-15,86255.2,-42.8148], Tmin=(934.877,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(743.318,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(174.604,'J/(mol*K)'), comment="""Thermo group additivity estimation: group(O2s-(Cds-Cd)H) + group(Cs-(Cds-Cds)CsHH) + group(Cs-CsHHH) + group(Cds-CdsCsH) + group(Cds-CdsOsH) + radical(CCJ2_triplet) + radical(Allyl_S) + radical(C=COJ) + radical(Cds_S)"""),
)

species(
    label = '[CH]=[C][CH]C1[C]=CO1(26185)',
    structure = SMILES('[CH]=[C][CH]C1[C]=CO1'),
    E0 = (835.648,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
    molecularWeight = (92.0954,'amu'),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[1.31841,0.0307041,6.51976e-05,-1.28152e-07,5.69903e-11,100628,24.5141], Tmin=(100,'K'), Tmax=(928.54,'K')), NASAPolynomial(coeffs=[26.5271,-0.00516459,5.65664e-06,-1.05237e-09,5.93527e-14,92811.3,-112.116], Tmin=(928.54,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(835.648,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(249.434,'J/(mol*K)'), comment="""Thermo group additivity estimation: group(O2s-Cs(Cds-Cd)) + group(Cs-(Cds-Cds)CsOsH) + group(Cs-(Cds-Cds)CsHH) + group(Cds-CdsCsH) + group(Cds-CdsCsH) + group(Cds-CdsOsH) + group(Cds-CdsHH) + ring(Cyclobutene) + radical(Cds_P) + radical(Cds_S) + radical(Cds_S) + radical(C=CCJC(O)C=C)"""),
)

species(
    label = 'C#CC=[C][C]=C[O](26186)',
    structure = SMILES('C#CC=[C][C]=C[O]'),
    E0 = (653.14,'kJ/mol'),
    modes = [
        HarmonicOscillator(frequencies=([1670,1700,300,440,750,770,3400,2100,2995,3025,975,1000,1300,1375,400,500,1630,1680,2175,525,180,180],'cm^-1')),
        HinderedRotor(inertia=(1.75351,'amu*angstrom^2'), symmetry=1, barrier=(40.3167,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(1.75332,'amu*angstrom^2'), symmetry=1, barrier=(40.3123,'kJ/mol'), semiclassical=False),
    ],
    spinMultiplicity = 4,
    opticalIsomers = 1,
    molecularWeight = (91.0874,'amu'),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[0.126342,0.0703187,-8.15493e-05,4.42774e-08,-8.84442e-12,78706.9,23.8862], Tmin=(100,'K'), Tmax=(1431.13,'K')), NASAPolynomial(coeffs=[19.5805,0.00349408,1.54057e-06,-5.07392e-10,4.08156e-14,74413.6,-72.4882], Tmin=(1431.13,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(653.14,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(224.491,'J/(mol*K)'), comment="""Thermo group additivity estimation: group(O2s-(Cds-Cd)H) + group(Cds-Cds(Cds-Cds)H) + group(Cds-Cds(Cds-Cds)H) + group(Cds-CdsCtH) + group(Cds-CdsOsH) + group(Ct-Ct(Cds-Cds)) + group(Ct-CtH) + radical(C=CJC=C) + radical(C=CJC=C) + radical(C=COJ)"""),
)

species(
    label = '[CH]=[C]C=[C]C=C[O](26187)',
    structure = SMILES('[CH]=[C]C=[C]C=C[O]'),
    E0 = (723.45,'kJ/mol'),
    modes = [
        HarmonicOscillator(frequencies=([1670,1700,300,440,2995,3010,3025,975,987.5,1000,1300,1337.5,1375,400,450,500,1630,1655,1680,3120,650,792.5,1650,180,180],'cm^-1')),
        HinderedRotor(inertia=(1.61857,'amu*angstrom^2'), symmetry=1, barrier=(37.214,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(1.62449,'amu*angstrom^2'), symmetry=1, barrier=(37.3502,'kJ/mol'), semiclassical=False),
    ],
    spinMultiplicity = 5,
    opticalIsomers = 1,
    molecularWeight = (92.0954,'amu'),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[-0.20668,0.0756015,-8.76164e-05,4.76389e-08,-9.48916e-12,87177,26.5889], Tmin=(100,'K'), Tmax=(1450.01,'K')), NASAPolynomial(coeffs=[20.494,0.00368415,2.10348e-06,-6.56437e-10,5.22281e-14,82730.9,-75.6021], Tmin=(1450.01,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(723.45,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(249.434,'J/(mol*K)'), comment="""Thermo group additivity estimation: group(O2s-(Cds-Cd)H) + group(Cds-Cds(Cds-Cds)H) + group(Cds-Cds(Cds-Cds)H) + group(Cds-Cds(Cds-Cds)H) + group(Cds-Cds(Cds-Cds)H) + group(Cds-CdsOsH) + group(Cds-CdsHH) + radical(Cds_P) + radical(C=CJC=C) + radical(C=COJ) + radical(C=CJC=C)"""),
)

species(
    label = '[CH]C=C=C[C]=C[O](26188)',
    structure = SMILES('[CH]C=C=C[C]=C[O]'),
    E0 = (700.824,'kJ/mol'),
    modes = [
        HarmonicOscillator(frequencies=([1685,370,2995,3010,3025,975,987.5,1000,1300,1337.5,1375,400,450,500,1630,1655,1680,540,610,2055,180,180,180,180,180],'cm^-1')),
        HinderedRotor(inertia=(2.02905,'amu*angstrom^2'), symmetry=1, barrier=(46.6518,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(2.02894,'amu*angstrom^2'), symmetry=1, barrier=(46.6494,'kJ/mol'), semiclassical=False),
    ],
    spinMultiplicity = 5,
    opticalIsomers = 1,
    molecularWeight = (92.0954,'amu'),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[0.0672002,0.0706144,-6.87667e-05,3.32427e-08,-6.11102e-12,84444.9,26.7593], Tmin=(100,'K'), Tmax=(1500.78,'K')), NASAPolynomial(coeffs=[18.5248,0.0131201,-3.00693e-06,3.46415e-10,-1.7325e-14,79839.4,-66.6673], Tmin=(1500.78,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(700.824,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(249.434,'J/(mol*K)'), comment="""Thermo group additivity estimation: group(O2s-(Cds-Cd)H) + group(Cs-(Cds-Cds)HHH) + group(Cds-CdsCsH) + group(Cds-Cds(Cds-Cds)H) + group(Cds-Cds(Cds-Cds)H) + group(Cds-CdsOsH) + group(Cdd-CdsCds) + radical(C=CJC=C) + radical(C=COJ) + radical(AllylJ2_triplet)"""),
)

species(
    label = '[CH]=CC=[C][C]=C[O](26189)',
    structure = SMILES('[CH]=CC=[C][C]=C[O]'),
    E0 = (723.45,'kJ/mol'),
    modes = [
        HarmonicOscillator(frequencies=([1670,1700,300,440,2995,3010,3025,975,987.5,1000,1300,1337.5,1375,400,450,500,1630,1655,1680,3120,650,792.5,1650,180,180],'cm^-1')),
        HinderedRotor(inertia=(1.61857,'amu*angstrom^2'), symmetry=1, barrier=(37.214,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(1.62449,'amu*angstrom^2'), symmetry=1, barrier=(37.3502,'kJ/mol'), semiclassical=False),
    ],
    spinMultiplicity = 5,
    opticalIsomers = 1,
    molecularWeight = (92.0954,'amu'),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[-0.206679,0.0756015,-8.76164e-05,4.76389e-08,-9.48915e-12,87177,26.5889], Tmin=(100,'K'), Tmax=(1450.01,'K')), NASAPolynomial(coeffs=[20.494,0.00368417,2.10347e-06,-6.56435e-10,5.2228e-14,82730.9,-75.602], Tmin=(1450.01,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(723.45,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(249.434,'J/(mol*K)'), comment="""Thermo group additivity estimation: group(O2s-(Cds-Cd)H) + group(Cds-Cds(Cds-Cds)H) + group(Cds-Cds(Cds-Cds)H) + group(Cds-Cds(Cds-Cds)H) + group(Cds-Cds(Cds-Cds)H) + group(Cds-CdsOsH) + group(Cds-CdsHH) + radical(C=COJ) + radical(C=CJC=C) + radical(Cds_P) + radical(C=CJC=C)"""),
)

species(
    label = 'C=[C]C=[C][C]=C[O](24989)',
    structure = SMILES('C=[C]C=[C][C]=C[O]'),
    E0 = (675.35,'kJ/mol'),
    modes = [
        HarmonicOscillator(frequencies=([2950,3100,1380,975,1025,1650,2995,3025,975,1000,1300,1375,400,500,1630,1680,1670,1685,1700,300,370,440,180,180,180],'cm^-1')),
        HinderedRotor(inertia=(1.77597,'amu*angstrom^2'), symmetry=1, barrier=(40.8331,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(1.78786,'amu*angstrom^2'), symmetry=1, barrier=(41.1063,'kJ/mol'), semiclassical=False),
    ],
    spinMultiplicity = 5,
    opticalIsomers = 1,
    molecularWeight = (92.0954,'amu'),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[-0.000283518,0.0762297,-9.20979e-05,5.33849e-08,-1.14246e-11,81380,25.0358], Tmin=(100,'K'), Tmax=(1336.68,'K')), NASAPolynomial(coeffs=[18.8718,0.00648097,1.06916e-06,-5.11427e-10,4.51732e-14,77520.7,-67.0517], Tmin=(1336.68,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(675.35,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(249.434,'J/(mol*K)'), comment="""Thermo group additivity estimation: group(O2s-(Cds-Cd)H) + group(Cds-Cds(Cds-Cds)H) + group(Cds-Cds(Cds-Cds)H) + group(Cds-Cds(Cds-Cds)H) + group(Cds-Cds(Cds-Cds)H) + group(Cds-CdsOsH) + group(Cds-CdsHH) + radical(C=CJC=C) + radical(C=CJC=C) + radical(C=COJ) + radical(C=CJC=C)"""),
)

species(
    label = '[CH]=C=C[C]=[C][CH]O(26190)',
    structure = SMILES('[CH]=C=C[C]=[C][CH]O'),
    E0 = (753.948,'kJ/mol'),
    modes = [
        HarmonicOscillator(frequencies=([3120,650,792.5,1650,3615,1277.5,1000,3025,407.5,1350,352.5,540,610,2055,1670,1700,300,440,3010,987.5,1337.5,450,1655,180],'cm^-1')),
        HinderedRotor(inertia=(1.58632,'amu*angstrom^2'), symmetry=1, barrier=(36.4726,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(1.58815,'amu*angstrom^2'), symmetry=1, barrier=(36.5146,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(1.59195,'amu*angstrom^2'), symmetry=1, barrier=(36.6022,'kJ/mol'), semiclassical=False),
    ],
    spinMultiplicity = 5,
    opticalIsomers = 1,
    molecularWeight = (92.0954,'amu'),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[0.927217,0.0640286,-6.73294e-05,3.00783e-08,-2.80772e-12,90793.1,26.743], Tmin=(100,'K'), Tmax=(861.021,'K')), NASAPolynomial(coeffs=[15.2888,0.0122104,-3.01483e-06,3.80325e-10,-2.07056e-14,87767.6,-43.6017], Tmin=(861.021,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(753.948,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(245.277,'J/(mol*K)'), comment="""Thermo group additivity estimation: group(O2s-CsH) + group(Cs-(Cds-Cds)OsHH) + group(Cds-CdsCsH) + group(Cds-Cds(Cds-Cds)H) + group(Cds-Cds(Cds-Cds)H) + group(Cds-CdsHH) + group(Cdd-CdsCds) + radical(Cds_S) + radical(C=CJC=C) + radical(C=CCJO) + radical(C=C=CJ)"""),
)

species(
    label = '[CH]C=C[CH][C]=C=O(26191)',
    structure = SMILES('[CH]C=C[CH][C]=C=O'),
    E0 = (640.336,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
    molecularWeight = (92.0954,'amu'),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[1.09558,0.0581781,-4.72783e-05,1.99749e-08,-3.42994e-12,77123.9,23.6575], Tmin=(100,'K'), Tmax=(1378.19,'K')), NASAPolynomial(coeffs=[13.0595,0.0234551,-9.48694e-06,1.69452e-09,-1.13975e-13,73826.2,-37.8993], Tmin=(1378.19,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(640.336,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(245.277,'J/(mol*K)'), comment="""Thermo group additivity estimation: group(Cs-Cd(CCO)HH) + group(Cs-(Cds-Cds)HHH) + group(Cds-CdsCsH) + group(Cds-CdsCsH) + group(Cds-(Cdd-O2d)CsH) + radical(C=CCJC=C=O) + radical(CCCJ=C=O) + radical(AllylJ2_triplet)"""),
)

species(
    label = '[CH]C=CC#CC=O(26192)',
    structure = SMILES('[CH]C=CC#CC=O'),
    E0 = (475.405,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
    molecularWeight = (92.0954,'amu'),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[1.22056,0.0526211,-3.62498e-05,1.21079e-08,-1.61733e-12,57285.3,24.634], Tmin=(100,'K'), Tmax=(1744.27,'K')), NASAPolynomial(coeffs=[15.2544,0.0204385,-8.57434e-06,1.53032e-09,-1.0129e-13,52389.5,-50.879], Tmin=(1744.27,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(475.405,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(245.277,'J/(mol*K)'), comment="""Thermo group additivity estimation: group(Cs-(Cds-Cds)HHH) + group(Cds-CdsCsH) + group(Cds-CdsCtH) + group(Cds-O2d(Cds-Cds)H) + group(Ct-Ct(Cds-Cds)) + group(Ct-CtCs) + radical(AllylJ2_triplet)"""),
)

species(
    label = '[CH]=C=C=CC=C[O](26193)',
    structure = SMILES('[CH]=C=C=CC=C[O]'),
    E0 = (462.225,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
    molecularWeight = (92.0954,'amu'),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[-0.218809,0.0696047,-7.35468e-05,3.62989e-08,-6.54161e-12,55764.7,26.7672], Tmin=(100,'K'), Tmax=(1624.8,'K')), NASAPolynomial(coeffs=[20.0588,0.00393998,1.60939e-06,-5.02105e-10,3.8405e-14,51253.6,-74.5073], Tmin=(1624.8,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(462.225,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(253.591,'J/(mol*K)'), comment="""Thermo group additivity estimation: group(O2s-(Cds-Cd)H) + group(Cds-Cds(Cds-Cds)H) + group(Cds-Cds(Cds-Cds)H) + group(Cds-CdsOsH) + group(Cds-CdsHH) + group(Cdd-CdsCds) + group(Cdd-CdsCds) + radical(C=COJ) + radical(C=C=CJ)"""),
)

species(
    label = 'C#CC=[C][C]=CO(26194)',
    structure = SMILES('C#CC=[C][C]=CO'),
    E0 = (511.677,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
    molecularWeight = (92.0954,'amu'),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[-0.534513,0.0802777,-9.56622e-05,5.23727e-08,-1.04076e-11,61720.6,24.8448], Tmin=(100,'K'), Tmax=(1468.86,'K')), NASAPolynomial(coeffs=[22.1485,0.00073289,3.72014e-06,-9.71836e-10,7.37827e-14,56974.4,-86.7822], Tmin=(1468.86,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(511.677,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(245.277,'J/(mol*K)'), comment="""Thermo group additivity estimation: group(O2s-(Cds-Cd)H) + group(Cds-Cds(Cds-Cds)H) + group(Cds-Cds(Cds-Cds)H) + group(Cds-CdsCtH) + group(Cds-CdsOsH) + group(Ct-Ct(Cds-Cds)) + group(Ct-CtH) + radical(C=CJC=C) + radical(C=CJC=C)"""),
)

species(
    label = '[CH]=C1C=C[C]1C=O(26082)',
    structure = SMILES('[CH]=C1C=C[C]1C=O'),
    E0 = (442.366,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
    molecularWeight = (92.0954,'amu'),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[1.85414,0.0304892,3.30491e-05,-7.08901e-08,3.08259e-11,53296.8,17.3828], Tmin=(100,'K'), Tmax=(951.996,'K')), NASAPolynomial(coeffs=[16.2348,0.0118121,-3.29958e-06,6.26734e-10,-5.11441e-14,48667,-61.2239], Tmin=(951.996,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(442.366,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(253.591,'J/(mol*K)'), comment="""Thermo group additivity estimation: group(Cs-(Cds-O2d)CsCsH) + group(Cds-Cds(Cds-Cds)Cs) + group(Cds-CdsCsH) + group(Cds-Cds(Cds-Cds)H) + group(Cds-OdCsH) + group(Cds-CdsHH) + ring(Cyclobutene) + radical(Cds_P) + radical(C=CCJ(C=O)C=C)"""),
)

species(
    label = '[CH]C=C=[CH](5318)',
    structure = SMILES('[CH]C=C=[CH]'),
    E0 = (671.864,'kJ/mol'),
    modes = [
        HarmonicOscillator(frequencies=([3120,650,792.5,1650,3010,987.5,1337.5,450,1655,540,610,2055,180,180],'cm^-1')),
        HinderedRotor(inertia=(2.15394,'amu*angstrom^2'), symmetry=1, barrier=(49.5233,'kJ/mol'), semiclassical=False),
    ],
    spinMultiplicity = 4,
    opticalIsomers = 1,
    molecularWeight = (51.0666,'amu'),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[2.71958,0.0243611,-3.13576e-06,-1.22602e-08,6.62132e-12,80856.1,14.3569], Tmin=(100,'K'), Tmax=(921.186,'K')), NASAPolynomial(coeffs=[6.85187,0.0161268,-5.53695e-06,9.18944e-10,-6.04195e-14,79682.8,-7.47571], Tmin=(921.186,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(671.864,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(153.818,'J/(mol*K)'), comment="""Thermo group additivity estimation: group(Cs-(Cds-Cds)HHH) + group(Cds-CdsCsH) + group(Cds-CdsHH) + group(Cdd-CdsCds) + radical(C=C=CJ) + radical(AllylJ2_triplet)"""),
)

species(
    label = '[C]C=O(24998)',
    structure = SMILES('[C]C=O'),
    E0 = (491.572,'kJ/mol'),
    modes = [
        HarmonicOscillator(frequencies=([2782.5,750,1395,475,1775,1000],'cm^-1')),
    ],
    spinMultiplicity = 4,
    opticalIsomers = 1,
    molecularWeight = (41.0287,'amu'),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[3.66873,0.00225797,1.94571e-05,-2.81427e-08,1.12145e-11,59138.9,7.7037], Tmin=(100,'K'), Tmax=(946.546,'K')), NASAPolynomial(coeffs=[6.76418,0.00196444,-3.42231e-07,7.48207e-11,-7.93789e-15,57980.1,-10.086], Tmin=(946.546,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(491.572,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(83.1447,'J/(mol*K)'), comment="""Thermo group additivity estimation: group(Cs-(Cds-O2d)HHH) + group(Cds-OdCsH) + radical(CJ3)"""),
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
    label = '[CH]C=C=C[O](24997)',
    structure = SMILES('[CH]C=C=C[O]'),
    E0 = (450.057,'kJ/mol'),
    modes = [
        HarmonicOscillator(frequencies=([2995,3025,975,1000,1300,1375,400,500,1630,1680,540,610,2055,180,180,180,180],'cm^-1')),
        HinderedRotor(inertia=(2.07021,'amu*angstrom^2'), symmetry=1, barrier=(47.5981,'kJ/mol'), semiclassical=False),
    ],
    spinMultiplicity = 4,
    opticalIsomers = 1,
    molecularWeight = (67.066,'amu'),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[2.08058,0.0321668,3.59383e-06,-3.27891e-08,1.63889e-11,54207.6,17.2931], Tmin=(100,'K'), Tmax=(940.194,'K')), NASAPolynomial(coeffs=[12.873,0.0113851,-3.3498e-06,5.6759e-10,-4.11617e-14,51067.4,-40.0161], Tmin=(940.194,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(450.057,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(178.761,'J/(mol*K)'), comment="""Thermo group additivity estimation: group(O2s-(Cds-Cd)H) + group(Cs-(Cds-Cds)HHH) + group(Cds-CdsCsH) + group(Cds-CdsOsH) + group(Cdd-CdsCds) + radical(AllylJ2_triplet) + radical(C=COJ)"""),
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
    label = '[CH]=C=C[CH][C]1[CH]O1(26195)',
    structure = SMILES('[CH]=C=C[CH][C]1[CH]O1'),
    E0 = (752.59,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
    molecularWeight = (92.0954,'amu'),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[0.315022,0.0728065,-8.94086e-05,5.475e-08,-1.24059e-11,90655.6,24.933], Tmin=(100,'K'), Tmax=(1280.8,'K')), NASAPolynomial(coeffs=[15.6897,0.0106415,-3.40897e-08,-3.95064e-10,4.14057e-14,87877.8,-48.5155], Tmin=(1280.8,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(752.59,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(249.434,'J/(mol*K)'), comment="""Thermo group additivity estimation: group(O2s-CsCs) + group(Cs-CsCsOsH) + group(Cs-(Cds-Cds)CsHH) + group(Cs-CsOsHH) + group(Cds-CdsCsH) + group(Cds-CdsHH) + group(Cdd-CdsCds) + ring(Ethylene_oxide) + radical(C=C=CJ) + radical(CCsJO) + radical(C=CCJCO) + radical(C2CsJO)"""),
)

species(
    label = '[CH]=C1[CH]C=[C]C1[O](26119)',
    structure = SMILES('[CH]=C1[CH]C=[C]C1[O]'),
    E0 = (774.226,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
    molecularWeight = (92.0954,'amu'),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[2.29146,0.0210987,5.00264e-05,-8.14822e-08,3.27938e-11,93194.3,22.8743], Tmin=(100,'K'), Tmax=(964.335,'K')), NASAPolynomial(coeffs=[13.5891,0.0154323,-5.23841e-06,1.02276e-09,-7.96904e-14,89099.8,-41.1516], Tmin=(964.335,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(774.226,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(257.749,'J/(mol*K)'), comment="""Thermo group additivity estimation: group(O2s-CsH) + group(Cs-(Cds-Cds)(Cds-Cds)OsH) + group(Cs-(Cds-Cds)(Cds-Cds)HH) + group(Cds-CdsCsCs) + group(Cds-CdsCsH) + group(Cds-CdsCsH) + group(Cds-CdsHH) + ring(4-Methylenecyclopentene) + radical(C=CCJC=C) + radical(Cds_P) + radical(CC(C)OJ) + radical(cyclopentene-vinyl)"""),
)

species(
    label = '[CH]C=CC=C=C=O(26196)',
    structure = SMILES('[CH]C=CC=C=C=O'),
    E0 = (545.925,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
    molecularWeight = (92.0954,'amu'),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[2.05822,0.0331537,3.87228e-06,-2.97074e-08,1.40204e-11,65738.1,7.21077], Tmin=(100,'K'), Tmax=(972.925,'K')), NASAPolynomial(coeffs=[11.6014,0.0167108,-5.9168e-06,1.07872e-09,-7.74559e-14,62802.4,-44.1116], Tmin=(972.925,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(545.925,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(249.434,'J/(mol*K)'), comment="""Thermo group additivity estimation: group(Cs-(Cds-Cds)HHH) + group(Cds-CdsCsH) + group(Cds-Cds(Cds-Cds)H) + group(Cds-Cds(Cds-Cds)H) + group(Cdd-CdsCds) + radical(AllylJ2_triplet)"""),
)

species(
    label = '[CH][C]=C[C]=C[CH][O](26197)',
    structure = SMILES('[CH][C]=C[C]=C[CH][O]'),
    E0 = (1019.26,'kJ/mol'),
    modes = [
        HarmonicOscillator(frequencies=([3025,407.5,1350,352.5,2995,3025,975,1000,1300,1375,400,500,1630,1680,1670,1700,300,440,378.548,378.572,378.59,378.592,378.601,378.617],'cm^-1')),
        HinderedRotor(inertia=(0.506059,'amu*angstrom^2'), symmetry=1, barrier=(51.4717,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.506103,'amu*angstrom^2'), symmetry=1, barrier=(51.472,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.506029,'amu*angstrom^2'), symmetry=1, barrier=(51.4716,'kJ/mol'), semiclassical=False),
    ],
    spinMultiplicity = 7,
    opticalIsomers = 1,
    molecularWeight = (92.0954,'amu'),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[1.4404,0.0581233,-5.35392e-05,2.88105e-08,-6.65218e-12,122679,27.1938], Tmin=(100,'K'), Tmax=(1015.04,'K')), NASAPolynomial(coeffs=[8.31802,0.0310196,-1.34849e-05,2.50253e-09,-1.72446e-13,121283,-6.08916], Tmin=(1015.04,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(1019.26,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(245.277,'J/(mol*K)'), comment="""Thermo group additivity estimation: group(O2s-CsH) + group(Cs-(Cds-Cds)OsHH) + group(Cs-(Cds-Cds)HHH) + group(Cds-CdsCsH) + group(Cds-CdsCsH) + group(Cds-Cds(Cds-Cds)H) + group(Cds-Cds(Cds-Cds)H) + radical(CCOJ) + radical(C=CCJO) + radical(Cds_S) + radical(C=CJC=C) + radical(AllylJ2_triplet)"""),
)

species(
    label = '[CH][C]=[C]C[C]=C[O](26198)',
    structure = SMILES('[CH][C]=[C]C[C]=C[O]'),
    E0 = (1069.6,'kJ/mol'),
    modes = [
        HarmonicOscillator(frequencies=([2750,2850,1437.5,1250,1305,750,350,3010,987.5,1337.5,450,1655,1670,1685,1700,300,370,440,278.694,278.736,278.818,278.822,278.839,278.846],'cm^-1')),
        HinderedRotor(inertia=(0.909765,'amu*angstrom^2'), symmetry=1, barrier=(50.1647,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.910111,'amu*angstrom^2'), symmetry=1, barrier=(50.1634,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.908978,'amu*angstrom^2'), symmetry=1, barrier=(50.1614,'kJ/mol'), semiclassical=False),
    ],
    spinMultiplicity = 7,
    opticalIsomers = 1,
    molecularWeight = (92.0954,'amu'),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[1.02196,0.0642828,-6.51317e-05,3.52114e-08,-7.69632e-12,128752,29.12], Tmin=(100,'K'), Tmax=(1103.06,'K')), NASAPolynomial(coeffs=[12.2343,0.0236238,-9.8416e-06,1.79528e-09,-1.22823e-13,126278,-26.0728], Tmin=(1103.06,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(1069.6,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(245.277,'J/(mol*K)'), comment="""Thermo group additivity estimation: group(O2s-(Cds-Cd)H) + group(Cs-(Cds-Cds)(Cds-Cds)HH) + group(Cs-(Cds-Cds)HHH) + group(Cds-CdsCsH) + group(Cds-CdsCsH) + group(Cds-CdsCsH) + group(Cds-CdsOsH) + radical(Cds_S) + radical(Cds_S) + radical(C=COJ) + radical(AllylJ2_triplet) + radical(Cds_S)"""),
)

species(
    label = '[CH][C]=C[CH]C=[C][O](26199)',
    structure = SMILES('[CH][C]=C[CH]C=[C][O]'),
    E0 = (935.39,'kJ/mol'),
    modes = [
        HarmonicOscillator(frequencies=([3025,407.5,1350,352.5,2995,3025,975,1000,1300,1375,400,500,1630,1680,1670,1700,300,440,461.897,461.897,461.897,461.898,461.898,461.899],'cm^-1')),
        HinderedRotor(inertia=(0.338771,'amu*angstrom^2'), symmetry=1, barrier=(51.2891,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.338771,'amu*angstrom^2'), symmetry=1, barrier=(51.289,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.338771,'amu*angstrom^2'), symmetry=1, barrier=(51.2891,'kJ/mol'), semiclassical=False),
    ],
    spinMultiplicity = 7,
    opticalIsomers = 1,
    molecularWeight = (92.0954,'amu'),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[1.45711,0.0463698,-1.17129e-05,-1.99764e-08,1.17097e-11,112602,29.4072], Tmin=(100,'K'), Tmax=(963.696,'K')), NASAPolynomial(coeffs=[12.9307,0.0212148,-7.53067e-06,1.32313e-09,-9.17493e-14,109347,-30.9348], Tmin=(963.696,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(935.39,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(245.277,'J/(mol*K)'), comment="""Thermo group additivity estimation: group(O2s-(Cds-Cd)H) + group(Cs-(Cds-Cds)(Cds-Cds)HH) + group(Cs-(Cds-Cds)HHH) + group(Cds-CdsCsH) + group(Cds-CdsCsH) + group(Cds-CdsCsH) + group(Cds-CdsOsH) + radical(AllylJ2_triplet) + radical(Cds_S) + radical(C=CCJC=C) + radical(C=CJO) + radical(C=COJ)"""),
)

species(
    label = '[CH][C]=[C][CH]C=C[O](26200)',
    structure = SMILES('[CH][C]=[C][CH]C=C[O]'),
    E0 = (933.487,'kJ/mol'),
    modes = [
        HarmonicOscillator(frequencies=([3025,407.5,1350,352.5,2995,3025,975,1000,1300,1375,400,500,1630,1680,1670,1700,300,440,353.601,354.652,356.557,356.851,357.606,357.738],'cm^-1')),
        HinderedRotor(inertia=(0.547568,'amu*angstrom^2'), symmetry=1, barrier=(50.2606,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.564556,'amu*angstrom^2'), symmetry=1, barrier=(50.1485,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.566853,'amu*angstrom^2'), symmetry=1, barrier=(50.2085,'kJ/mol'), semiclassical=False),
    ],
    spinMultiplicity = 7,
    opticalIsomers = 1,
    molecularWeight = (92.0954,'amu'),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[1.22759,0.049056,-9.08839e-06,-3.01429e-08,1.72216e-11,112383,27.5762], Tmin=(100,'K'), Tmax=(938.024,'K')), NASAPolynomial(coeffs=[15.7912,0.0167411,-5.0486e-06,8.40919e-10,-5.90476e-14,108341,-48.7387], Tmin=(938.024,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(933.487,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(245.277,'J/(mol*K)'), comment="""Thermo group additivity estimation: group(O2s-(Cds-Cd)H) + group(Cs-(Cds-Cds)(Cds-Cds)HH) + group(Cs-(Cds-Cds)HHH) + group(Cds-CdsCsH) + group(Cds-CdsCsH) + group(Cds-CdsCsH) + group(Cds-CdsOsH) + radical(Cds_S) + radical(Cds_S) + radical(C=COJ) + radical(AllylJ2_triplet) + radical(C=CCJC=C)"""),
)

species(
    label = '[CH][C]=[C]C=[C]C[O](26201)',
    structure = SMILES('[CH][C]=[C]C=[C]C[O]'),
    E0 = (1139.8,'kJ/mol'),
    modes = [
        HarmonicOscillator(frequencies=([2750,2850,1437.5,1250,1305,750,350,3010,987.5,1337.5,450,1655,1670,1685,1700,300,370,440,379.13,379.139,379.148,379.155,379.159,379.183],'cm^-1')),
        HinderedRotor(inertia=(0.525653,'amu*angstrom^2'), symmetry=1, barrier=(53.6255,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.525682,'amu*angstrom^2'), symmetry=1, barrier=(53.6271,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.525741,'amu*angstrom^2'), symmetry=1, barrier=(53.6252,'kJ/mol'), semiclassical=False),
    ],
    spinMultiplicity = 7,
    opticalIsomers = 1,
    molecularWeight = (92.0954,'amu'),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[1.07111,0.0754159,-0.000124688,1.18057e-07,-4.25255e-11,137182,28.388], Tmin=(100,'K'), Tmax=(877.713,'K')), NASAPolynomial(coeffs=[3.32186,0.0391191,-1.81569e-05,3.34069e-09,-2.23122e-13,137789,23.5368], Tmin=(877.713,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(1139.8,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(245.277,'J/(mol*K)'), comment="""Thermo group additivity estimation: group(O2s-CsH) + group(Cs-(Cds-Cds)OsHH) + group(Cs-(Cds-Cds)HHH) + group(Cds-CdsCsH) + group(Cds-CdsCsH) + group(Cds-Cds(Cds-Cds)H) + group(Cds-Cds(Cds-Cds)H) + radical(CCOJ) + radical(C=CJC=C) + radical(Cds_S) + radical(AllylJ2_triplet) + radical(Cds_S)"""),
)

species(
    label = '[CH]=[C]C[C]=C=C[O](24793)',
    structure = SMILES('[CH]=[C]C[C]=C=C[O]'),
    E0 = (884.778,'kJ/mol'),
    modes = [
        HarmonicOscillator(frequencies=([3120,650,792.5,1650,2750,2850,1437.5,1250,1305,750,350,540,610,2055,3010,987.5,1337.5,450,1655,1670,1700,300,440,180,180],'cm^-1')),
        HinderedRotor(inertia=(1.02673,'amu*angstrom^2'), symmetry=1, barrier=(23.6066,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(1.02283,'amu*angstrom^2'), symmetry=1, barrier=(23.5169,'kJ/mol'), semiclassical=False),
    ],
    spinMultiplicity = 5,
    opticalIsomers = 1,
    molecularWeight = (92.0954,'amu'),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[0.874011,0.0646744,-7.26473e-05,4.05177e-08,-8.77503e-12,106530,27.5235], Tmin=(100,'K'), Tmax=(1136.36,'K')), NASAPolynomial(coeffs=[15.4021,0.0135351,-5.1432e-06,9.15066e-10,-6.23981e-14,103228,-44.4233], Tmin=(1136.36,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(884.778,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(249.434,'J/(mol*K)'), comment="""Thermo group additivity estimation: group(O2s-(Cds-Cd)H) + group(Cs-(Cds-Cds)(Cds-Cds)HH) + group(Cds-CdsCsH) + group(Cds-CdsCsH) + group(Cds-CdsHH) + group(Cds-CdsOsH) + group(Cdd-CdsCds) + radical(Cds_S) + radical(C=COJ) + radical(Cds_S) + radical(Cds_P)"""),
)

species(
    label = '[CH]=[C]C1C=[C]C1[O](26202)',
    structure = SMILES('[CH]=[C]C1C=[C]C1[O]'),
    E0 = (1011.15,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
    molecularWeight = (92.0954,'amu'),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[1.32004,0.0525694,-4.65932e-05,2.06555e-08,-3.61998e-12,121715,24.2966], Tmin=(100,'K'), Tmax=(1376.82,'K')), NASAPolynomial(coeffs=[14.0434,0.0156053,-6.32222e-06,1.15608e-09,-7.93485e-14,118211,-41.1547], Tmin=(1376.82,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(1011.15,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(253.591,'J/(mol*K)'), comment="""Thermo group additivity estimation: group(O2s-CsH) + group(Cs-(Cds-Cds)(Cds-Cds)CsH) + group(Cs-(Cds-Cds)CsOsH) + group(Cds-CdsCsH) + group(Cds-CdsCsH) + group(Cds-CdsCsH) + group(Cds-CdsHH) + ring(Cyclobutene) + radical(cyclobutene-vinyl) + radical(Cds_S) + radical(CC(C)OJ) + radical(Cds_P)"""),
)

species(
    label = '[O]C1[C]=CC=[C][CH]1(26203)',
    structure = SMILES('[O]C1[C]=CC=[C][CH]1'),
    E0 = (691.754,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
    molecularWeight = (92.0954,'amu'),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[2.15062,0.0185214,7.16541e-05,-1.12698e-07,4.5741e-11,83285.5,20.729], Tmin=(100,'K'), Tmax=(953.805,'K')), NASAPolynomial(coeffs=[17.6008,0.00965886,-2.36879e-06,5.21372e-10,-4.91865e-14,77794,-66.4155], Tmin=(953.805,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(691.754,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(257.749,'J/(mol*K)'), comment="""Thermo group additivity estimation: group(O2s-CsH) + group(Cs-(Cds-Cds)CsOsH) + group(Cs-(Cds-Cds)CsHH) + group(Cds-CdsCsH) + group(Cds-CdsCsH) + group(Cds-Cds(Cds-Cds)H) + group(Cds-Cds(Cds-Cds)H) + ring(1,3-Cyclohexadiene) + radical(Cds_S) + radical(Cds_S) + radical(C=CCJC(O)C=C) + radical(CC(C)OJ)"""),
)

species(
    label = '[CH]=[C]CC=C=C=O(24786)',
    structure = SMILES('[CH]=[C]CC=C=C=O'),
    E0 = (691.033,'kJ/mol'),
    modes = [
        HarmonicOscillator(frequencies=([3120,650,792.5,1650,540,610,2055,2750,2850,1437.5,1250,1305,750,350,2120,512.5,787.5,1685,370,3010,987.5,1337.5,450,1655,4000],'cm^-1')),
        HinderedRotor(inertia=(0.000525832,'amu*angstrom^2'), symmetry=1, barrier=(5.97033,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.000524487,'amu*angstrom^2'), symmetry=1, barrier=(5.95505,'kJ/mol'), semiclassical=False),
    ],
    spinMultiplicity = 3,
    opticalIsomers = 1,
    molecularWeight = (92.0954,'amu'),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[2.18904,0.0364105,-2.87925e-05,1.16452e-08,-1.89495e-12,83180,9.53559], Tmin=(100,'K'), Tmax=(1453.2,'K')), NASAPolynomial(coeffs=[10.3176,0.014036,-5.69729e-06,1.05001e-09,-7.22e-14,80817.5,-32.7183], Tmin=(1453.2,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(691.033,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(249.434,'J/(mol*K)'), comment="""Thermo group additivity estimation: group(Cs-(Cds-Cds)(Cds-Cds)HH) + group(Cds-CdsCsH) + group(Cds-CdsCsH) + group(Cds-CdsHH) + group(Cdd-CdsCds) + radical(Cds_S) + radical(Cds_P)"""),
)

species(
    label = '[CH]=C=C[C]=[C]C[O](26204)',
    structure = SMILES('[CH]=C=C[C]=[C]C[O]'),
    E0 = (862.358,'kJ/mol'),
    modes = [
        HarmonicOscillator(frequencies=([3120,650,792.5,1650,2750,2850,1437.5,1250,1305,750,350,540,610,2055,3010,987.5,1337.5,450,1655,1670,1700,300,440,180,180],'cm^-1')),
        HinderedRotor(inertia=(2.65711,'amu*angstrom^2'), symmetry=1, barrier=(61.0922,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.650782,'amu*angstrom^2'), symmetry=1, barrier=(14.9628,'kJ/mol'), semiclassical=False),
    ],
    spinMultiplicity = 5,
    opticalIsomers = 1,
    molecularWeight = (92.0954,'amu'),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[1.1005,0.0721754,-0.000119368,1.06011e-07,-3.56405e-11,103814,26.3807], Tmin=(100,'K'), Tmax=(907.117,'K')), NASAPolynomial(coeffs=[6.55591,0.02763,-1.18276e-05,2.07622e-09,-1.33789e-13,103667,5.23983], Tmin=(907.117,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(862.358,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(249.434,'J/(mol*K)'), comment="""Thermo group additivity estimation: group(O2s-CsH) + group(Cs-(Cds-Cds)OsHH) + group(Cds-CdsCsH) + group(Cds-Cds(Cds-Cds)H) + group(Cds-Cds(Cds-Cds)H) + group(Cds-CdsHH) + group(Cdd-CdsCds) + radical(CCOJ) + radical(C=CJC=C) + radical(Cds_S) + radical(C=C=CJ)"""),
)

species(
    label = '[CH]=C=[C]C=[C]C[O](26205)',
    structure = SMILES('[CH]=C=[C]C=[C]C[O]'),
    E0 = (901.204,'kJ/mol'),
    modes = [
        HarmonicOscillator(frequencies=([3120,650,792.5,1650,2750,2850,1437.5,1250,1305,750,350,540,610,2055,3010,987.5,1337.5,450,1655,1670,1700,300,440,180,180],'cm^-1')),
        HinderedRotor(inertia=(0.536135,'amu*angstrom^2'), symmetry=1, barrier=(12.3268,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.537464,'amu*angstrom^2'), symmetry=1, barrier=(12.3573,'kJ/mol'), semiclassical=False),
    ],
    spinMultiplicity = 5,
    opticalIsomers = 1,
    molecularWeight = (92.0954,'amu'),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[1.17182,0.0699288,-0.000114064,1.01254e-07,-3.44406e-11,108484,27.0015], Tmin=(100,'K'), Tmax=(884.436,'K')), NASAPolynomial(coeffs=[6.73371,0.0273376,-1.22569e-05,2.22317e-09,-1.47075e-13,108182,4.70702], Tmin=(884.436,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(901.204,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(249.434,'J/(mol*K)'), comment="""Thermo group additivity estimation: group(O2s-CsH) + group(Cs-(Cds-Cds)OsHH) + group(Cds-CdsCsH) + group(Cds-Cds(Cds-Cds)H) + group(Cds-Cds(Cds-Cds)H) + group(Cds-CdsHH) + group(Cdd-CdsCds) + radical(Cds_S) + radical(C=C=CJ) + radical(CCOJ) + radical(Cds_S)"""),
)

species(
    label = '[CH]=[C]C[C][C]=C[O](26206)',
    structure = SMILES('[CH]=[C]C[C][C]=C[O]'),
    E0 = (1167.36,'kJ/mol'),
    modes = [
        HarmonicOscillator(frequencies=([3120,650,792.5,1650,2750,2850,1437.5,1250,1305,750,350,3010,987.5,1337.5,450,1655,1670,1700,300,440,223.944,223.984,224.002,224.047],'cm^-1')),
        HinderedRotor(inertia=(0.603878,'amu*angstrom^2'), symmetry=1, barrier=(21.486,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.603145,'amu*angstrom^2'), symmetry=1, barrier=(21.4874,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.603759,'amu*angstrom^2'), symmetry=1, barrier=(21.4866,'kJ/mol'), semiclassical=False),
    ],
    spinMultiplicity = 7,
    opticalIsomers = 1,
    molecularWeight = (92.0954,'amu'),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[0.730937,0.0680341,-7.92342e-05,4.54174e-08,-1.00697e-11,140523,29.3002], Tmin=(100,'K'), Tmax=(1112.61,'K')), NASAPolynomial(coeffs=[16.0707,0.0128846,-4.8818e-06,8.6554e-10,-5.89267e-14,137109,-46.3423], Tmin=(1112.61,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(1167.36,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(245.277,'J/(mol*K)'), comment="""Thermo group additivity estimation: group(O2s-(Cds-Cd)H) + group(Cs-(Cds-Cds)CsHH) + group(Cs-(Cds-Cds)CsHH) + group(Cds-CdsCsH) + group(Cds-CdsCsH) + group(Cds-CdsOsH) + group(Cds-CdsHH) + radical(Cds_S) + radical(C=COJ) + radical(CCJ2_triplet) + radical(Cds_P) + radical(Cds_S)"""),
)

species(
    label = 'C#CC=C[C]=C[O](22635)',
    structure = SMILES('C#CC=C[C]=C[O]'),
    E0 = (454.145,'kJ/mol'),
    modes = [
        HarmonicOscillator(frequencies=([1685,370,750,770,3400,2100,2995,3010,3025,975,987.5,1000,1300,1337.5,1375,400,450,500,1630,1655,1680,2175,525,180,180],'cm^-1')),
        HinderedRotor(inertia=(1.62614,'amu*angstrom^2'), symmetry=1, barrier=(37.3881,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(1.62728,'amu*angstrom^2'), symmetry=1, barrier=(37.4144,'kJ/mol'), semiclassical=False),
    ],
    spinMultiplicity = 3,
    opticalIsomers = 1,
    molecularWeight = (92.0954,'amu'),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[0.831519,0.0562972,-2.5861e-05,-2.38426e-08,1.8337e-11,54747.6,21.7645], Tmin=(100,'K'), Tmax=(911.449,'K')), NASAPolynomial(coeffs=[20.957,0.00392216,1.17357e-06,-3.44429e-10,2.26331e-14,49585.7,-81.6544], Tmin=(911.449,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(454.145,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(249.434,'J/(mol*K)'), comment="""Thermo group additivity estimation: group(O2s-(Cds-Cd)H) + group(Cds-Cds(Cds-Cds)H) + group(Cds-Cds(Cds-Cds)H) + group(Cds-CdsCtH) + group(Cds-CdsOsH) + group(Ct-Ct(Cds-Cds)) + group(Ct-CtH) + radical(C=CJC=C) + radical(C=COJ)"""),
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
    label = '[C][CH]C=C=[CH](21432)',
    structure = SMILES('[C][CH]C=C=[CH]'),
    E0 = (1089.12,'kJ/mol'),
    modes = [
        HarmonicOscillator(frequencies=([3120,650,792.5,1650,3025,407.5,1350,352.5,3010,987.5,1337.5,450,1655,540,610,2055,180],'cm^-1')),
        HinderedRotor(inertia=(1.89172,'amu*angstrom^2'), symmetry=1, barrier=(43.4943,'kJ/mol'), semiclassical=False),
    ],
    spinMultiplicity = 6,
    opticalIsomers = 1,
    molecularWeight = (63.0773,'amu'),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[2.21402,0.0332244,-1.71911e-05,-7.36501e-09,7.12796e-12,131061,15.3217], Tmin=(100,'K'), Tmax=(920.935,'K')), NASAPolynomial(coeffs=[11.6622,0.00802887,-1.95572e-06,2.84609e-10,-1.9169e-14,128649,-33.1293], Tmin=(920.935,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(1089.12,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(178.761,'J/(mol*K)'), comment="""Thermo group additivity estimation: group(Cs-(Cds-Cds)CsHH) + group(Cs-CsHHH) + group(Cds-CdsCsH) + group(Cds-CdsHH) + group(Cdd-CdsCds) + radical(CJ3) + radical(C=C=CJ) + radical(Allyl_S)"""),
)

species(
    label = 'C#CC=[C]C=C[O](26207)',
    structure = SMILES('C#CC=[C]C=C[O]'),
    E0 = (454.145,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
    molecularWeight = (92.0954,'amu'),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[0.831519,0.0562972,-2.5861e-05,-2.38426e-08,1.8337e-11,54747.6,21.7645], Tmin=(100,'K'), Tmax=(911.449,'K')), NASAPolynomial(coeffs=[20.957,0.00392216,1.17357e-06,-3.44429e-10,2.26331e-14,49585.7,-81.6544], Tmin=(911.449,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(454.145,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(249.434,'J/(mol*K)'), comment="""Thermo group additivity estimation: group(O2s-(Cds-Cd)H) + group(Cds-Cds(Cds-Cds)H) + group(Cds-Cds(Cds-Cds)H) + group(Cds-CdsCtH) + group(Cds-CdsOsH) + group(Ct-Ct(Cds-Cds)) + group(Ct-CtH) + radical(C=CJC=C) + radical(C=COJ)"""),
)

species(
    label = 'C#CC1[CH]C1=C[O](26162)',
    structure = SMILES('C#CC1[CH]C1=C[O]'),
    E0 = (473.979,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
    molecularWeight = (92.0954,'amu'),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[1.17742,0.0465123,-1.46167e-06,-4.5491e-08,2.47975e-11,57122.6,17.4421], Tmin=(100,'K'), Tmax=(927.36,'K')), NASAPolynomial(coeffs=[19.9989,0.00539895,2.26244e-07,-1.11479e-10,3.34677e-15,51908.7,-81.2312], Tmin=(927.36,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(473.979,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(253.591,'J/(mol*K)'), comment="""Thermo group additivity estimation: group(O2s-(Cds-Cd)H) + group(Cs-(Cds-Cds)CtCsH) + group(Cs-(Cds-Cds)CsHH) + group(Cds-CdsCsCs) + group(Cds-CdsOsH) + group(Ct-CtCs) + group(Ct-CtH) + ring(Methylene_cyclopropane) + radical(Allyl_S) + radical(C=COJ)"""),
)

species(
    label = 'C#C[CH]C1[C]C1[O](26208)',
    structure = SMILES('C#C[CH]C1[C]C1[O]'),
    E0 = (907.892,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
    molecularWeight = (92.0954,'amu'),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[1.35083,0.047119,-1.29485e-05,-2.96095e-08,1.9115e-11,109300,22.4113], Tmin=(100,'K'), Tmax=(894.63,'K')), NASAPolynomial(coeffs=[16.9558,0.00846893,-3.25849e-07,-1.3137e-10,1.143e-14,105262,-58.0967], Tmin=(894.63,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(907.892,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(249.434,'J/(mol*K)'), comment="""Thermo group additivity estimation: group(O2s-CsH) + group(Cs-CsCsCsH) + group(Cs-CsCsOsH) + group(Cs-CsCsHH) + group(Cs-CtCsHH) + group(Ct-CtCs) + group(Ct-CtH) + ring(Cyclopropane) + radical(CC(C)OJ) + radical(Sec_Propargyl) + radical(CCJ2_triplet)"""),
)

species(
    label = '[C]#CC=C[C]C[O](26209)',
    structure = SMILES('[C]#CC=C[C]C[O]'),
    E0 = (1063.27,'kJ/mol'),
    modes = [
        HarmonicOscillator(frequencies=([2750,2850,1437.5,1250,1305,750,350,2995,3025,975,1000,1300,1375,400,500,1630,1680,2175,525,250.014,250.26,250.355,250.436,1316.55],'cm^-1')),
        HinderedRotor(inertia=(0.256334,'amu*angstrom^2'), symmetry=1, barrier=(11.3688,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(1.40991,'amu*angstrom^2'), symmetry=1, barrier=(62.5825,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(1.41041,'amu*angstrom^2'), symmetry=1, barrier=(62.5876,'kJ/mol'), semiclassical=False),
    ],
    spinMultiplicity = 5,
    opticalIsomers = 1,
    molecularWeight = (92.0954,'amu'),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[1.20524,0.066807,-9.82816e-05,8.27937e-08,-2.79595e-11,127978,25.8514], Tmin=(100,'K'), Tmax=(825.516,'K')), NASAPolynomial(coeffs=[7.69394,0.0275908,-1.28956e-05,2.42826e-09,-1.6623e-13,127171,-2.60398], Tmin=(825.516,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(1063.27,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(245.277,'J/(mol*K)'), comment="""Thermo group additivity estimation: group(O2s-CsH) + group(Cs-(Cds-Cds)CsHH) + group(Cs-CsOsHH) + group(Cds-CdsCsH) + group(Cds-CdsCtH) + group(Ct-Ct(Cds-Cds)) + group(Ct-CtH) + radical(CCOJ) + radical(Acetyl) + radical(CCJ2_triplet)"""),
)

species(
    label = '[C]#CC[CH][C][CH][O](26210)',
    structure = SMILES('[C]#CC[CH][C][CH][O]'),
    E0 = (1328.01,'kJ/mol'),
    modes = [
        HarmonicOscillator(frequencies=([2750,2850,1437.5,1250,1305,750,350,2175,525,3000,3050,390,425,1340,1360,335,370,200,800,1000,1200,1400,1600],'cm^-1')),
        HinderedRotor(inertia=(0.156089,'amu*angstrom^2'), symmetry=1, barrier=(3.5888,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.156089,'amu*angstrom^2'), symmetry=1, barrier=(3.5888,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.156089,'amu*angstrom^2'), symmetry=1, barrier=(3.5888,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.156089,'amu*angstrom^2'), symmetry=1, barrier=(3.5888,'kJ/mol'), semiclassical=False),
    ],
    spinMultiplicity = 7,
    opticalIsomers = 1,
    molecularWeight = (92.0954,'amu'),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[0.917934,0.0862496,-0.000176117,1.7595e-07,-6.31478e-11,159817,29.7975], Tmin=(100,'K'), Tmax=(910.361,'K')), NASAPolynomial(coeffs=[1.67906,0.0362437,-1.68387e-05,3.00519e-09,-1.92718e-13,161612,36.8167], Tmin=(910.361,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(1328.01,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(241.12,'J/(mol*K)'), comment="""Thermo group additivity estimation: group(O2s-CsH) + group(Cs-CsCsHH) + group(Cs-CsCsHH) + group(Cs-CtCsHH) + group(Cs-CsOsHH) + group(Ct-CtCs) + group(Ct-CtH) + radical(CCJ2_triplet) + radical(Acetyl) + radical(RCCJCC) + radical(CCsJOH) + radical(CCOJ)"""),
)

species(
    label = '[C]#C[CH]C[C][CH][O](26211)',
    structure = SMILES('[C]#C[CH]C[C][CH][O]'),
    E0 = (1279.77,'kJ/mol'),
    modes = [
        HarmonicOscillator(frequencies=([2750,2850,1437.5,1250,1305,750,350,2175,525,3000,3050,390,425,1340,1360,335,370,200,800,1000,1200,1400,1600],'cm^-1')),
        HinderedRotor(inertia=(0.156089,'amu*angstrom^2'), symmetry=1, barrier=(3.5888,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.156089,'amu*angstrom^2'), symmetry=1, barrier=(3.5888,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.156089,'amu*angstrom^2'), symmetry=1, barrier=(3.5888,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.156089,'amu*angstrom^2'), symmetry=1, barrier=(3.5888,'kJ/mol'), semiclassical=False),
    ],
    spinMultiplicity = 7,
    opticalIsomers = 1,
    molecularWeight = (92.0954,'amu'),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[0.650333,0.0867719,-0.000161276,1.48424e-07,-5.00352e-11,154029,28.0322], Tmin=(100,'K'), Tmax=(930.506,'K')), NASAPolynomial(coeffs=[6.37902,0.0280878,-1.17733e-05,1.97616e-09,-1.20847e-13,154437,8.72999], Tmin=(930.506,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(1279.77,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(241.12,'J/(mol*K)'), comment="""Thermo group additivity estimation: group(O2s-CsH) + group(Cs-CsCsHH) + group(Cs-CsCsHH) + group(Cs-CtCsHH) + group(Cs-CsOsHH) + group(Ct-CtCs) + group(Ct-CtH) + radical(CCsJOH) + radical(Sec_Propargyl) + radical(CCJ2_triplet) + radical(CCOJ) + radical(Acetyl)"""),
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
    E0 = (656.042,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS2',
    E0 = (817.963,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS3',
    E0 = (759.958,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS4',
    E0 = (880.799,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS5',
    E0 = (820.312,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS6',
    E0 = (949.379,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS7',
    E0 = (969.845,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS8',
    E0 = (912.551,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS9',
    E0 = (837.547,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS10',
    E0 = (878.367,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS11',
    E0 = (916.883,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS12',
    E0 = (859.478,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS13',
    E0 = (847.288,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS14',
    E0 = (894.875,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS15',
    E0 = (1105.69,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS16',
    E0 = (1089.98,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS17',
    E0 = (854.178,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS18',
    E0 = (862.439,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS19',
    E0 = (785.106,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS20',
    E0 = (719.259,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS21',
    E0 = (678.903,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS22',
    E0 = (719.442,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS23',
    E0 = (658.971,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS24',
    E0 = (1073.32,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS25',
    E0 = (1368.36,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS26',
    E0 = (1134.25,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS27',
    E0 = (1236.59,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS28',
    E0 = (1142.12,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS29',
    E0 = (857.368,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS30',
    E0 = (1036.14,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS31',
    E0 = (816.491,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS32',
    E0 = (976.873,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS33',
    E0 = (832.716,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS34',
    E0 = (714.618,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS35',
    E0 = (786.75,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS36',
    E0 = (1329.57,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS37',
    E0 = (835.648,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS38',
    E0 = (881.626,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS39',
    E0 = (937.914,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS40',
    E0 = (863.603,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS41',
    E0 = (916.829,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS42',
    E0 = (955.395,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS43',
    E0 = (904.02,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS44',
    E0 = (689.082,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS45',
    E0 = (734.289,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS46',
    E0 = (745.01,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS47',
    E0 = (681.015,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS48',
    E0 = (664.326,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS49',
    E0 = (1197.75,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS50',
    E0 = (1331.72,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS51',
    E0 = (816.095,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS52',
    E0 = (1131.62,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS53',
    E0 = (843.241,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS54',
    E0 = (774.226,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS55',
    E0 = (681.015,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS56',
    E0 = (1042.12,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS57',
    E0 = (1097.46,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS58',
    E0 = (963.242,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS59',
    E0 = (996.887,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS60',
    E0 = (1164.78,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS61',
    E0 = (1062.46,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS62',
    E0 = (1011.49,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS63',
    E0 = (691.754,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS64',
    E0 = (691.033,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS65',
    E0 = (1055.44,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS66',
    E0 = (1051.28,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS67',
    E0 = (1190.23,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS68',
    E0 = (656.042,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS69',
    E0 = (962.515,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS70',
    E0 = (1156.8,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS71',
    E0 = (678.903,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS72',
    E0 = (658.971,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS73',
    E0 = (908.066,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS74',
    E0 = (1137.72,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS75',
    E0 = (1391.41,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS76',
    E0 = (1304.74,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

reaction(
    label = 'reaction1',
    reactants = ['[CH]=C=C[CH][C]=C[O](22638)'],
    products = ['C#CC=O(21959)', '[CH]=C=[CH](18734)'],
    transitionState = 'TS1',
    kinetics = Arrhenius(A=(5e+12,'s^-1'), n=0, Ea=(0,'kJ/mol'), T0=(1,'K'), Tmin=(300,'K'), Tmax=(1500,'K'), comment="""Exact match found for rate rule [RJJ]
Euclidian distance = 0
family: 1,4_Linear_birad_scission"""),
)

reaction(
    label = 'reaction2',
    reactants = ['[CH]=C=C[CH][C]=C[O](22638)'],
    products = ['[CH]=[C]C1[CH]C1=C[O](26150)'],
    transitionState = 'TS2',
    kinetics = Arrhenius(A=(1.10037e+11,'s^-1'), n=0.4025, Ea=(161.921,'kJ/mol'), T0=(1,'K'), comment="""Estimated using average of templates [R4;multiplebond_intra;radadd_intra_cddouble] + [R4;doublebond_intra;radadd_intra] for rate rule [R4;doublebond_intra;radadd_intra_cddouble]
Euclidian distance = 1.0
family: Intra_R_Add_Exocyclic"""),
)

reaction(
    label = 'reaction3',
    reactants = ['[CH]=C=C[CH][C]=C[O](22638)'],
    products = ['[CH]=[C]C1[CH][C]=CO1(25321)'],
    transitionState = 'TS3',
    kinetics = Arrhenius(A=(1.92333e+10,'s^-1'), n=0.385799, Ea=(103.916,'kJ/mol'), T0=(1,'K'), comment="""Estimated using average of templates [Rn;doublebond_intra;radadd_intra_O] + [R6;doublebond_intra;radadd_intra] for rate rule [R6;doublebond_intra;radadd_intra_O]
Euclidian distance = 1.0
family: Intra_R_Add_Exocyclic"""),
)

reaction(
    label = 'reaction4',
    reactants = ['H(8)', '[CH]=C=C=C[C]=C[O](26164)'],
    products = ['[CH]=C=C[CH][C]=C[O](22638)'],
    transitionState = 'TS4',
    kinetics = Arrhenius(A=(8.22e+08,'cm^3/(mol*s)'), n=1.533, Ea=(7.77387,'kJ/mol'), T0=(1,'K'), comment="""From training reaction 192 used for Cd_R;HJ
Exact match found for rate rule [Cd_R;HJ]
Euclidian distance = 0
family: R_Addition_MultipleBond"""),
)

reaction(
    label = 'reaction5',
    reactants = ['H(8)', '[CH]=C=C[CH][C]=C=O(26165)'],
    products = ['[CH]=C=C[CH][C]=C[O](22638)'],
    transitionState = 'TS5',
    kinetics = Arrhenius(A=(8.22e+08,'cm^3/(mol*s)'), n=1.533, Ea=(7.77387,'kJ/mol'), T0=(1,'K'), comment="""From training reaction 192 used for Cd_R;HJ
Exact match found for rate rule [Cd_R;HJ]
Euclidian distance = 0
family: R_Addition_MultipleBond"""),
)

reaction(
    label = 'reaction6',
    reactants = ['O(T)(63)', '[CH]=C=CC=C=[CH](21250)'],
    products = ['[CH]=C=C[CH][C]=C[O](22638)'],
    transitionState = 'TS6',
    kinetics = Arrhenius(A=(4.61703e+06,'m^3/(mol*s)'), n=-0.199588, Ea=(22.3125,'kJ/mol'), T0=(1,'K'), comment="""Estimated using average of templates [R_R;O_atom_triplet] + [Ct_Ct;YJ] for rate rule [Ct_Ct;O_atom_triplet]
Euclidian distance = 2.0
Multiplied by reaction path degeneracy 2.0
family: R_Addition_MultipleBond"""),
)

reaction(
    label = 'reaction7',
    reactants = ['[CH]=C=[C]C[C]=C[O](26166)'],
    products = ['[CH]=C=C[CH][C]=C[O](22638)'],
    transitionState = 'TS7',
    kinetics = Arrhenius(A=(7.80481e+09,'s^-1'), n=1.09833, Ea=(177.686,'kJ/mol'), T0=(1,'K'), comment="""Estimated using an average for rate rule [R2H_S;Cd_rad_out_double;XH_out]
Euclidian distance = 0
Multiplied by reaction path degeneracy 2.0
family: intra_H_migration"""),
)

reaction(
    label = 'reaction8',
    reactants = ['[CH]=C=C[CH][C]=[C]O(26167)'],
    products = ['[CH]=C=C[CH][C]=C[O](22638)'],
    transitionState = 'TS8',
    kinetics = Arrhenius(A=(6.25466e+06,'s^-1'), n=1.80084, Ea=(158.227,'kJ/mol'), T0=(1,'K'), comment="""Estimated using average of templates [R2H_S;Y_rad_out;O_H_out] + [R2H_S;Cd_rad_out;XH_out] for rate rule [R2H_S;Cd_rad_out;O_H_out]
Euclidian distance = 1.0
family: intra_H_migration"""),
)

reaction(
    label = 'reaction9',
    reactants = ['[CH]=C=C[CH][C]=C[O](22638)'],
    products = ['[CH]=C=C[CH]C=[C][O](26168)'],
    transitionState = 'TS9',
    kinetics = Arrhenius(A=(2.97782e+08,'s^-1'), n=1.48417, Ea=(181.505,'kJ/mol'), T0=(1,'K'), comment="""Estimated using an average for rate rule [R2H_D;Cd_rad_out;XH_out]
Euclidian distance = 0
family: intra_H_migration"""),
)

reaction(
    label = 'reaction10',
    reactants = ['[CH]=C=[C][CH]C=C[O](26169)'],
    products = ['[CH]=C=C[CH][C]=C[O](22638)'],
    transitionState = 'TS10',
    kinetics = Arrhenius(A=(191.5,'s^-1'), n=3.05, Ea=(222.325,'kJ/mol'), T0=(1,'K'), comment="""Estimated using template [R3Hall;Cd_rad_out;Cd_H_out_doubleC] for rate rule [R3HJ;Cd_rad_out_double;Cd_H_out_doubleC]
Euclidian distance = 1.41421356237
family: intra_H_migration"""),
)

reaction(
    label = 'reaction11',
    reactants = ['C#C[CH]C[C][C]=O(26170)'],
    products = ['[CH]=C=C[CH][C]=C[O](22638)'],
    transitionState = 'TS11',
    kinetics = Arrhenius(A=(1.34586e+06,'s^-1'), n=1.99734, Ea=(143.275,'kJ/mol'), T0=(1,'K'), comment="""Estimated using template [R3Hall;Y_rad_out;Cs_H_out_H/Cd] for rate rule [R3HJ;Cd_rad_out;Cs_H_out_H/Cd]
Euclidian distance = 1.41421356237
Multiplied by reaction path degeneracy 2.0
family: intra_H_migration"""),
)

reaction(
    label = 'reaction12',
    reactants = ['C=[C][C]=C[C]=C[O](24988)'],
    products = ['[CH]=C=C[CH][C]=C[O](22638)'],
    transitionState = 'TS12',
    kinetics = Arrhenius(A=(8.97864e+06,'s^-1'), n=1.84533, Ea=(184.129,'kJ/mol'), T0=(1,'K'), comment="""Estimated using average of templates [R3H;Y_rad_out;Cd_H_out_singleH] + [R3H;Cd_rad_out;Cd_H_out_single] for rate rule [R3H;Cd_rad_out;Cd_H_out_singleH]
Euclidian distance = 1.0
Multiplied by reaction path degeneracy 2.0
family: intra_H_migration"""),
)

reaction(
    label = 'reaction13',
    reactants = ['[CH]=C=[C][CH][C]=CO(26171)'],
    products = ['[CH]=C=C[CH][C]=C[O](22638)'],
    transitionState = 'TS13',
    kinetics = Arrhenius(A=(6.00274e+06,'s^-1'), n=1.61277, Ea=(94.8673,'kJ/mol'), T0=(1,'K'), comment="""Estimated using average of templates [RnH;Cd_rad_out_double;XH_out] + [R5Hall;Cd_rad_out;XH_out] for rate rule [R5Hall;Cd_rad_out_double;O_H_out]
Euclidian distance = 1.41421356237
family: intra_H_migration"""),
)

reaction(
    label = 'reaction14',
    reactants = ['[CH2][C]=C[CH][C]=C=O(24990)'],
    products = ['[CH]=C=C[CH][C]=C[O](22638)'],
    transitionState = 'TS14',
    kinetics = Arrhenius(A=(2.36537e+10,'s^-1'), n=1.06641, Ea=(235.882,'kJ/mol'), T0=(1,'K'), comment="""Estimated using template [RnH;Cd_rad_out;Cd_H_out_singleH] for rate rule [R6Hall;Cd_rad_out;Cd_H_out_singleH]
Euclidian distance = 1.0
Multiplied by reaction path degeneracy 2.0
family: intra_H_migration"""),
)

reaction(
    label = 'reaction15',
    reactants = ['H(8)', '[CH]=C=[C][CH][C]=C[O](26172)'],
    products = ['[CH]=C=C[CH][C]=C[O](22638)'],
    transitionState = 'TS15',
    kinetics = Arrhenius(A=(4.34078e+06,'m^3/(mol*s)'), n=0.278577, Ea=(0,'kJ/mol'), T0=(1,'K'), comment="""Estimated using an average for rate rule [Y_rad;H_rad]
Euclidian distance = 0
family: R_Recombination
Ea raised from -1.4 to 0 kJ/mol."""),
)

reaction(
    label = 'reaction16',
    reactants = ['H(8)', '[CH][C]=C[CH][C]=C=O(26173)'],
    products = ['[CH]=C=C[CH][C]=C[O](22638)'],
    transitionState = 'TS16',
    kinetics = Arrhenius(A=(4.34078e+06,'m^3/(mol*s)'), n=0.278577, Ea=(0,'kJ/mol'), T0=(1,'K'), comment="""Estimated using an average for rate rule [Y_rad;H_rad]
Euclidian distance = 0
family: R_Recombination
Ea raised from -1.4 to 0 kJ/mol."""),
)

reaction(
    label = 'reaction17',
    reactants = ['[CH]=C=C[CH][C]=C[O](22638)'],
    products = ['[CH]=C1[CH]C1[C]=C[O](26085)'],
    transitionState = 'TS17',
    kinetics = Arrhenius(A=(4.00692e+11,'s^-1'), n=0.347401, Ea=(198.136,'kJ/mol'), T0=(1,'K'), comment="""Estimated using an average for rate rule [R3_D;doublebond_intra;radadd_intra_cs]
Euclidian distance = 0
family: Intra_R_Add_Endocyclic"""),
)

reaction(
    label = 'reaction18',
    reactants = ['[CH]=C=C[CH][C]=C[O](22638)'],
    products = ['[O]C=[C][CH]C1[C]=C1(26174)'],
    transitionState = 'TS18',
    kinetics = Arrhenius(A=(4.2354e+12,'s^-1'), n=-0.1205, Ea=(206.398,'kJ/mol'), T0=(1,'K'), comment="""Estimated using average of templates [Rn;doublebond_intra_CdCdd;radadd_intra_cdsingleH] + [R3;doublebond_intra_CdCdd;radadd_intra_cdsingle] for rate rule [R3;doublebond_intra_CdCdd;radadd_intra_cdsingleH]
Euclidian distance = 1.0
family: Intra_R_Add_Endocyclic"""),
)

reaction(
    label = 'reaction19',
    reactants = ['[CH]=C=C[CH][C]=C[O](22638)'],
    products = ['[CH]C1=C[CH]C1=C[O](26175)'],
    transitionState = 'TS19',
    kinetics = Arrhenius(A=(2.77275e+09,'s^-1'), n=0.708429, Ea=(129.065,'kJ/mol'), T0=(1,'K'), comment="""Estimated using template [R4_linear;doublebond_intra;radadd_intra] for rate rule [R4_linear;doublebond_intra;radadd_intra_cddouble]
Euclidian distance = 1.0
family: Intra_R_Add_Endocyclic"""),
)

reaction(
    label = 'reaction20',
    reactants = ['[CH]=C=C[CH][C]=C[O](22638)'],
    products = ['[CH]C1=CC=[C][CH]O1(26110)'],
    transitionState = 'TS20',
    kinetics = Arrhenius(A=(2.49515e+10,'s^-1'), n=0.243684, Ea=(63.2167,'kJ/mol'), T0=(1,'K'), comment="""Estimated using template [R6_linear;doublebond_intra;radadd_intra] for rate rule [R6_linear;doublebond_intra;radadd_intra_O]
Euclidian distance = 1.0
family: Intra_R_Add_Endocyclic"""),
)

reaction(
    label = 'reaction21',
    reactants = ['[CH]=C=C[CH][C]=C[O](22638)'],
    products = ['[CH]=C=C[CH]C=C=O(26176)'],
    transitionState = 'TS21',
    kinetics = Arrhenius(A=(1.949e+11,'s^-1'), n=0.486, Ea=(22.8614,'kJ/mol'), T0=(1,'K'), comment="""Estimated using an average for rate rule [R2radExo;Y_rad;XH_Rrad]
Euclidian distance = 0
family: Intra_Disproportionation"""),
)

reaction(
    label = 'reaction22',
    reactants = ['[CH]=C=C[CH][C]=C[O](22638)'],
    products = ['[CH2]C#CC=C=C[O](24976)'],
    transitionState = 'TS22',
    kinetics = Arrhenius(A=(7.437e+08,'s^-1'), n=1.045, Ea=(63.4002,'kJ/mol'), T0=(1,'K'), comment="""Estimated using an average for rate rule [R3radExo;Y_rad;XH_Rrad]
Euclidian distance = 0
family: Intra_Disproportionation"""),
)

reaction(
    label = 'reaction23',
    reactants = ['[CH]=C=C[CH][C]=C[O](22638)'],
    products = ['[CH]=C=CC=C1[CH]O1(26177)'],
    transitionState = 'TS23',
    kinetics = Arrhenius(A=(2.18e+16,'s^-1'), n=0, Ea=(2.9288,'kJ/mol'), T0=(1,'K'), comment="""Estimated using template [R3;Y_rad_out;Ypri_rad_out] for rate rule [R3_SD;O_rad;Ypri_rad_out]
Euclidian distance = 1.41421356237
family: Birad_recombination"""),
)

reaction(
    label = 'reaction24',
    reactants = ['[CH][C]=C[O](21209)', '[CH]=C=[CH](18734)'],
    products = ['[CH]=C=C[CH][C]=C[O](22638)'],
    transitionState = 'TS24',
    kinetics = Arrhenius(A=(2.29708e+06,'m^3/(mol*s)'), n=0.575199, Ea=(34.3157,'kJ/mol'), T0=(1,'K'), comment="""Estimated using template [Y_rad;Birad] for rate rule [Cd_pri_rad;Birad]
Euclidian distance = 2.0
Multiplied by reaction path degeneracy 2.0
family: Birad_R_Recombination"""),
)

reaction(
    label = 'reaction25',
    reactants = ['O(T)(63)', '[CH][C]=CC=C=[CH](21252)'],
    products = ['[CH]=C=C[CH][C]=C[O](22638)'],
    transitionState = 'TS25',
    kinetics = Arrhenius(A=(187219,'m^3/(mol*s)'), n=1.13083, Ea=(163.847,'kJ/mol'), T0=(1,'K'), Tmin=(303.03,'K'), Tmax=(2000,'K'), comment="""From training reaction 2 used for Y_rad;O_birad
Exact match found for rate rule [Y_rad;O_birad]
Euclidian distance = 0
Multiplied by reaction path degeneracy 2.0
family: Birad_R_Recombination"""),
)

reaction(
    label = 'reaction26',
    reactants = ['H(8)', '[CH]=[C]C=[C][C]=C[O](26178)'],
    products = ['[CH]=C=C[CH][C]=C[O](22638)'],
    transitionState = 'TS26',
    kinetics = Arrhenius(A=(1e+07,'m^3/(mol*s)'), n=0, Ea=(0,'kJ/mol'), T0=(1,'K'), comment="""Estimated using an average for rate rule [H_rad;Birad]
Euclidian distance = 0
family: Birad_R_Recombination"""),
)

reaction(
    label = 'reaction27',
    reactants = ['H(8)', '[C]#C[CH][CH][C]=C[O](26179)'],
    products = ['[CH]=C=C[CH][C]=C[O](22638)'],
    transitionState = 'TS27',
    kinetics = Arrhenius(A=(1e+07,'m^3/(mol*s)'), n=0, Ea=(0,'kJ/mol'), T0=(1,'K'), comment="""Estimated using an average for rate rule [H_rad;Birad]
Euclidian distance = 0
family: Birad_R_Recombination"""),
)

reaction(
    label = 'reaction28',
    reactants = ['[C]=[C]CC=C=C[O](24806)'],
    products = ['[CH]=C=C[CH][C]=C[O](22638)'],
    transitionState = 'TS28',
    kinetics = Arrhenius(A=(2.30972e+07,'s^-1'), n=1.70216, Ea=(184.173,'kJ/mol'), T0=(1,'K'), comment="""Estimated using template [R3H_MS;Y_rad_out;XH_out] for rate rule [R3H_TS;Ct_rad_out;XH_out]
Euclidian distance = 1.41421356237
Multiplied by reaction path degeneracy 2.0
family: intra_H_migration"""),
)

reaction(
    label = 'reaction29',
    reactants = ['[CH]=[C]C[CH][C]=C=O(24794)'],
    products = ['[CH]=C=C[CH][C]=C[O](22638)'],
    transitionState = 'TS29',
    kinetics = Arrhenius(A=(74200,'s^-1'), n=2.23, Ea=(44.3086,'kJ/mol'), T0=(1,'K'), comment="""Estimated using template [R4Hall;Cd_rad_out;Cs_H_out_H/OneDe] for rate rule [R4Hall;Cd_rad_out;Cs_H_out_H/Ct]
Euclidian distance = 1.0
Multiplied by reaction path degeneracy 2.0
family: intra_H_migration"""),
)

reaction(
    label = 'reaction30',
    reactants = ['[C]#C[CH]C[C]=C[O](26180)'],
    products = ['[CH]=C=C[CH][C]=C[O](22638)'],
    transitionState = 'TS30',
    kinetics = Arrhenius(A=(645397,'s^-1'), n=1.96429, Ea=(152.474,'kJ/mol'), T0=(1,'K'), comment="""Estimated using template [R4Hall;Y_rad_out;XH_out] for rate rule [R4HJ_2;Ct_rad_out;XH_out]
Euclidian distance = 1.41421356237
Multiplied by reaction path degeneracy 2.0
family: intra_H_migration"""),
)

reaction(
    label = 'reaction31',
    reactants = ['[C]#C[CH][CH]C=C[O](26181)'],
    products = ['[CH]=C=C[CH][C]=C[O](22638)'],
    transitionState = 'TS31',
    kinetics = Arrhenius(A=(1.658e+09,'s^-1'), n=0.699, Ea=(29.5516,'kJ/mol'), T0=(1,'K'), comment="""Estimated using template [R5Hall;Y_rad_out;Cd_H_out_doubleC] for rate rule [R5Hall;Ct_rad_out;Cd_H_out_doubleC]
Euclidian distance = 1.0
family: intra_H_migration"""),
)

reaction(
    label = 'reaction32',
    reactants = ['[C]#C[CH][CH][C]=CO(26182)'],
    products = ['[CH]=C=C[CH][C]=C[O](22638)'],
    transitionState = 'TS32',
    kinetics = Arrhenius(A=(1494.33,'s^-1'), n=2.17569, Ea=(93.5537,'kJ/mol'), T0=(1,'K'), comment="""Estimated using average of templates [RnH;Y_rad_out;O_H_out] + [R7Hall;Y_rad_out;XH_out] for rate rule [R7Hall;Ct_rad_out;O_H_out]
Euclidian distance = 1.41421356237
family: intra_H_migration"""),
)

reaction(
    label = 'reaction33',
    reactants = ['[CH]=C=C[CH][C]=C[O](22638)'],
    products = ['[O]C=[C]C1[CH][C]=C1(26183)'],
    transitionState = 'TS33',
    kinetics = Arrhenius(A=(9.88848e+07,'s^-1'), n=1.14903, Ea=(176.674,'kJ/mol'), T0=(1,'K'), comment="""Estimated using average of templates [Rn;triplebond_intra_H;radadd_intra_cs] + [R4_linear;multiplebond_intra;radadd_intra_cs] for rate rule [R4_linear;triplebond_intra_H;radadd_intra_cs]
Euclidian distance = 2.0
family: Intra_R_Add_Endocyclic"""),
)

reaction(
    label = 'reaction34',
    reactants = ['[CH]=C=C[CH][C]=C[O](22638)'],
    products = ['[O]C=C1[CH][C]=C[CH]1(26160)'],
    transitionState = 'TS34',
    kinetics = Arrhenius(A=(3.47e+11,'s^-1'), n=0.15, Ea=(58.576,'kJ/mol'), T0=(1,'K'), comment="""Estimated using an average for rate rule [R5_linear;triplebond_intra_H;radadd_intra_cddouble]
Euclidian distance = 0
family: Intra_R_Add_Endocyclic"""),
)

reaction(
    label = 'reaction35',
    reactants = ['[CH]=C=C[CH][C]=C[O](22638)'],
    products = ['[C]1[CH]OC=[C][CH]C=1(25334)'],
    transitionState = 'TS35',
    kinetics = Arrhenius(A=(2.07e+10,'s^-1'), n=0.124, Ea=(130.708,'kJ/mol'), T0=(1,'K'), comment="""Estimated using template [R7_linear;triplebond_intra_H;radadd_intra] for rate rule [R7_linear;triplebond_intra_H;radadd_intra_O]
Euclidian distance = 1.0
family: Intra_R_Add_Endocyclic"""),
)

reaction(
    label = 'reaction36',
    reactants = ['[C]#C(5143)', '[CH][CH][C]=C[O](26184)'],
    products = ['[CH]=C=C[CH][C]=C[O](22638)'],
    transitionState = 'TS36',
    kinetics = Arrhenius(A=(1.14854e+06,'m^3/(mol*s)'), n=0.575199, Ea=(34.3157,'kJ/mol'), T0=(1,'K'), comment="""Estimated using template [Y_rad;Birad] for rate rule [Ct_rad;Birad]
Euclidian distance = 1.0
family: Birad_R_Recombination"""),
)

reaction(
    label = 'reaction37',
    reactants = ['[CH]=C=C[CH][C]=C[O](22638)'],
    products = ['[CH]=[C][CH]C1[C]=CO1(26185)'],
    transitionState = 'TS37',
    kinetics = Arrhenius(A=(2.724e+10,'s^-1'), n=0.478, Ea=(179.606,'kJ/mol'), T0=(1,'K'), comment="""Estimated using an average for rate rule [R5;doublebond_intra;radadd_intra_O]
Euclidian distance = 0
family: Intra_R_Add_Exocyclic
Ea raised from 178.0 to 179.6 kJ/mol to match endothermicity of reaction."""),
)

reaction(
    label = 'reaction38',
    reactants = ['H(8)', 'C#CC=[C][C]=C[O](26186)'],
    products = ['[CH]=C=C[CH][C]=C[O](22638)'],
    transitionState = 'TS38',
    kinetics = Arrhenius(A=(6.23e+08,'cm^3/(mol*s)'), n=1.429, Ea=(16.6816,'kJ/mol'), T0=(1,'K'), comment="""From training reaction 193 used for Ca_Ca;HJ
Exact match found for rate rule [Ca_Ca;HJ]
Euclidian distance = 0
family: R_Addition_MultipleBond"""),
)

reaction(
    label = 'reaction39',
    reactants = ['[CH]=[C]C=[C]C=C[O](26187)'],
    products = ['[CH]=C=C[CH][C]=C[O](22638)'],
    transitionState = 'TS39',
    kinetics = Arrhenius(A=(1.806e+09,'s^-1'), n=1.172, Ea=(214.463,'kJ/mol'), T0=(1,'K'), comment="""Estimated using an average for rate rule [R2H_S;Cd_rad_out_Cd;Cd_H_out_doubleC]
Euclidian distance = 0
family: intra_H_migration"""),
)

reaction(
    label = 'reaction40',
    reactants = ['[CH]C=C=C[C]=C[O](26188)'],
    products = ['[CH]=C=C[CH][C]=C[O](22638)'],
    transitionState = 'TS40',
    kinetics = Arrhenius(A=(4.96519e+09,'s^-1'), n=1.05826, Ea=(162.779,'kJ/mol'), T0=(1,'K'), comment="""Estimated using an average for rate rule [R2H_S;Cd_rad_out_Cd;XH_out]
Euclidian distance = 0
family: intra_H_migration"""),
)

reaction(
    label = 'reaction41',
    reactants = ['[CH]=CC=[C][C]=C[O](26189)'],
    products = ['[CH]=C=C[CH][C]=C[O](22638)'],
    transitionState = 'TS41',
    kinetics = Arrhenius(A=(35492.9,'s^-1'), n=2.38823, Ea=(193.378,'kJ/mol'), T0=(1,'K'), comment="""Estimated using an average for rate rule [R3H_DS;Cd_rad_out;XH_out]
Euclidian distance = 0
family: intra_H_migration"""),
)

reaction(
    label = 'reaction42',
    reactants = ['C=[C]C=[C][C]=C[O](24989)'],
    products = ['[CH]=C=C[CH][C]=C[O](22638)'],
    transitionState = 'TS42',
    kinetics = Arrhenius(A=(2.19427e+11,'s^-1'), n=0.834916, Ea=(280.046,'kJ/mol'), T0=(1,'K'), comment="""Estimated using template [R4Hall;Cd_rad_out;Cd_H_out_singleH] for rate rule [R4HJ_2;Cd_rad_out;Cd_H_out_singleH]
Euclidian distance = 1.0
Multiplied by reaction path degeneracy 2.0
family: intra_H_migration"""),
)

reaction(
    label = 'reaction43',
    reactants = ['[CH]=C=C[C]=[C][CH]O(26190)'],
    products = ['[CH]=C=C[CH][C]=C[O](22638)'],
    transitionState = 'TS43',
    kinetics = Arrhenius(A=(1.13341e+09,'s^-1'), n=1.04717, Ea=(150.072,'kJ/mol'), T0=(1,'K'), comment="""Estimated using template [R4Hall;Cd_rad_out_Cd;XH_out] for rate rule [R4HJ_1;Cd_rad_out_Cd;O_H_out]
Euclidian distance = 1.41421356237
family: intra_H_migration"""),
)

reaction(
    label = 'reaction44',
    reactants = ['[CH]=C=C[CH][C]=C[O](22638)'],
    products = ['[CH]C=C[CH][C]=C=O(26191)'],
    transitionState = 'TS44',
    kinetics = Arrhenius(A=(136000,'s^-1'), n=1.9199, Ea=(33.0402,'kJ/mol'), T0=(1,'K'), comment="""Estimated using template [R5Hall;Cd_rad_out;XH_out] for rate rule [R5HJ_3;Cd_rad_out;XH_out]
Euclidian distance = 1.0
family: intra_H_migration"""),
)

reaction(
    label = 'reaction45',
    reactants = ['[CH]=C=C[CH][C]=C[O](22638)'],
    products = ['[CH]C=CC#CC=O(26192)'],
    transitionState = 'TS45',
    kinetics = Arrhenius(A=(2.00399e+09,'s^-1'), n=0.37, Ea=(78.2471,'kJ/mol'), T0=(1,'K'), comment="""Estimated using average of templates [R3;Y_rad;XH_Rrad_De] + [R3radExo;Y_rad;XH_Rrad] for rate rule [R3radExo;Y_rad;XH_Rrad_De]
Euclidian distance = 1.0
family: Intra_Disproportionation"""),
)

reaction(
    label = 'reaction46',
    reactants = ['[CH]=C=C[CH][C]=C[O](22638)'],
    products = ['[CH]=C=C=CC=C[O](26193)'],
    transitionState = 'TS46',
    kinetics = Arrhenius(A=(2.6374e+09,'s^-1'), n=0.37, Ea=(88.9686,'kJ/mol'), T0=(1,'K'), comment="""Estimated using average of templates [R3;Y_rad_De;XH_Rrad] + [R3radExo;Y_rad;XH_Rrad] for rate rule [R3radExo;Y_rad_De;XH_Rrad]
Euclidian distance = 1.0
family: Intra_Disproportionation"""),
)

reaction(
    label = 'reaction47',
    reactants = ['[CH]=C=C[CH][C]=C[O](22638)'],
    products = ['C#CC=[C][C]=CO(26194)'],
    transitionState = 'TS47',
    kinetics = Arrhenius(A=(5.14222e+08,'s^-1'), n=0.311, Ea=(24.9733,'kJ/mol'), T0=(1,'K'), comment="""Estimated using template [R4;Y_rad;XH_Rrad] for rate rule [R4radEndo;Y_rad;XH_Rrad]
Euclidian distance = 1.0
family: Intra_Disproportionation"""),
)

reaction(
    label = 'reaction48',
    reactants = ['[CH]=C=C[CH][C]=C[O](22638)'],
    products = ['[CH]=C1C=C[C]1C=O(26082)'],
    transitionState = 'TS48',
    kinetics = Arrhenius(A=(1.62e+12,'s^-1'), n=-0.305, Ea=(8.28432,'kJ/mol'), T0=(1,'K'), comment="""Estimated using template [R4;Y_rad_out;Ypri_rad_out] for rate rule [R4_SDS;Y_rad_out;Ypri_rad_out]
Euclidian distance = 1.0
family: Birad_recombination"""),
)

reaction(
    label = 'reaction49',
    reactants = ['[CH]C=C=[CH](5318)', '[C]C=O(24998)'],
    products = ['[CH]=C=C[CH][C]=C[O](22638)'],
    transitionState = 'TS49',
    kinetics = Arrhenius(A=(1.14854e+06,'m^3/(mol*s)'), n=0.575199, Ea=(34.3157,'kJ/mol'), T0=(1,'K'), comment="""Estimated using template [Y_rad;Birad] for rate rule [Cd_pri_rad;Birad]
Euclidian distance = 2.0
family: Birad_R_Recombination"""),
)

reaction(
    label = 'reaction50',
    reactants = ['[C]=[CH](18830)', '[CH]C=C=C[O](24997)'],
    products = ['[CH]=C=C[CH][C]=C[O](22638)'],
    transitionState = 'TS50',
    kinetics = Arrhenius(A=(1.14854e+06,'m^3/(mol*s)'), n=0.575199, Ea=(34.3157,'kJ/mol'), T0=(1,'K'), comment="""Estimated using template [Y_rad;Birad] for rate rule [Cd_pri_rad;Birad]
Euclidian distance = 2.0
family: Birad_R_Recombination"""),
)

reaction(
    label = 'reaction51',
    reactants = ['[CH]=C=C[O](8556)', '[CH]=C=[CH](18734)'],
    products = ['[CH]=C=C[CH][C]=C[O](22638)'],
    transitionState = 'TS51',
    kinetics = Arrhenius(A=(1.49708,'m^3/(mol*s)'), n=1.927, Ea=(54.5384,'kJ/mol'), T0=(1,'K'), comment="""Estimated using an average for rate rule [Ct_Ct;CdsJ=Cdd]
Euclidian distance = 0
Multiplied by reaction path degeneracy 2.0
family: R_Addition_MultipleBond"""),
)

reaction(
    label = 'reaction52',
    reactants = ['[CH]=C=C[O](8556)', '[CH][C]=[CH](21256)'],
    products = ['[CH]=C=C[CH][C]=C[O](22638)'],
    transitionState = 'TS52',
    kinetics = Arrhenius(A=(81155.9,'m^3/(mol*s)'), n=0.702818, Ea=(0,'kJ/mol'), T0=(1,'K'), comment="""Estimated using an average for rate rule [Y_rad;Cd_allenic]
Euclidian distance = 0
Multiplied by reaction path degeneracy 2.0
family: R_Recombination
Ea raised from -11.4 to 0 kJ/mol."""),
)

reaction(
    label = 'reaction53',
    reactants = ['[CH]=C=C[CH][C]=C[O](22638)'],
    products = ['[CH]=C=C[CH][C]1[CH]O1(26195)'],
    transitionState = 'TS53',
    kinetics = Arrhenius(A=(6.59543e+11,'s^-1'), n=0.400725, Ea=(187.199,'kJ/mol'), T0=(1,'K'), comment="""Estimated using template [R3_D;doublebond_intra;radadd_intra] for rate rule [R3_D;doublebond_intra;radadd_intra_O]
Euclidian distance = 1.0
family: Intra_R_Add_Endocyclic"""),
)

reaction(
    label = 'reaction54',
    reactants = ['[CH]=C=C[CH][C]=C[O](22638)'],
    products = ['[CH]=C1[CH]C=[C]C1[O](26119)'],
    transitionState = 'TS54',
    kinetics = Arrhenius(A=(6.01304e+12,'s^-1'), n=-0.3725, Ea=(118.184,'kJ/mol'), T0=(1,'K'), comment="""Estimated using an average for rate rule [R5_MS;doublebond_intra_CdCdd;radadd_intra]
Euclidian distance = 0
family: Intra_R_Add_Endocyclic
Ea raised from 114.7 to 118.2 kJ/mol to match endothermicity of reaction."""),
)

reaction(
    label = 'reaction55',
    reactants = ['[CH]=C=C[CH][C]=C[O](22638)'],
    products = ['[CH]C=CC=C=C=O(26196)'],
    transitionState = 'TS55',
    kinetics = Arrhenius(A=(2.1261e+09,'s^-1'), n=0.137, Ea=(24.9733,'kJ/mol'), T0=(1,'K'), comment="""Estimated using template [R5;Y_rad;XH_Rrad] for rate rule [R5radExo;Y_rad;XH_Rrad]
Euclidian distance = 1.0
family: Intra_Disproportionation"""),
)

reaction(
    label = 'reaction56',
    reactants = ['[CH][C]=C[C]=C[CH][O](26197)'],
    products = ['[CH]=C=C[CH][C]=C[O](22638)'],
    transitionState = 'TS56',
    kinetics = Arrhenius(A=(1.949e+11,'s^-1'), n=0.486, Ea=(22.8614,'kJ/mol'), T0=(1,'K'), comment="""Estimated using an average for rate rule [R2radExo;Y_rad_De;XH_Rrad]
Euclidian distance = 0
family: Intra_Disproportionation"""),
)

reaction(
    label = 'reaction57',
    reactants = ['[CH][C]=[C]C[C]=C[O](26198)'],
    products = ['[CH]=C=C[CH][C]=C[O](22638)'],
    transitionState = 'TS57',
    kinetics = Arrhenius(A=(2.94659e+10,'s^-1'), n=0.2847, Ea=(27.8529,'kJ/mol'), T0=(1,'K'), comment="""Estimated using average of templates [Rn;Y_rad;XH_Rrad_De] + [R2radExo;Y_rad;XH_Rrad] for rate rule [R2radExo;Y_rad;XH_Rrad_De]
Euclidian distance = 1.0
Multiplied by reaction path degeneracy 2.0
family: Intra_Disproportionation"""),
)

reaction(
    label = 'reaction58',
    reactants = ['[CH][C]=C[CH]C=[C][O](26199)'],
    products = ['[CH]=C=C[CH][C]=C[O](22638)'],
    transitionState = 'TS58',
    kinetics = Arrhenius(A=(1.4733e+10,'s^-1'), n=0.2847, Ea=(27.8529,'kJ/mol'), T0=(1,'K'), comment="""Estimated using average of templates [Rn;Y_rad;XH_Rrad_De] + [R2radExo;Y_rad;XH_Rrad] for rate rule [R2radExo;Y_rad;XH_Rrad_De]
Euclidian distance = 1.0
family: Intra_Disproportionation"""),
)

reaction(
    label = 'reaction59',
    reactants = ['[CH][C]=[C][CH]C=C[O](26200)'],
    products = ['[CH]=C=C[CH][C]=C[O](22638)'],
    transitionState = 'TS59',
    kinetics = Arrhenius(A=(7.437e+08,'s^-1'), n=1.045, Ea=(63.4002,'kJ/mol'), T0=(1,'K'), comment="""Estimated using an average for rate rule [R3radExo;Y_rad;XH_Rrad]
Euclidian distance = 0
family: Intra_Disproportionation"""),
)

reaction(
    label = 'reaction60',
    reactants = ['[CH][C]=[C]C=[C]C[O](26201)'],
    products = ['[CH]=C=C[CH][C]=C[O](22638)'],
    transitionState = 'TS60',
    kinetics = Arrhenius(A=(1.02844e+09,'s^-1'), n=0.311, Ea=(24.9733,'kJ/mol'), T0=(1,'K'), comment="""Estimated using template [R4;Y_rad;XH_Rrad] for rate rule [R4radEndo;Y_rad;XH_Rrad]
Euclidian distance = 1.0
Multiplied by reaction path degeneracy 2.0
family: Intra_Disproportionation"""),
)

reaction(
    label = 'reaction61',
    reactants = ['[CH]=[C]C[C]=C=C[O](24793)'],
    products = ['[CH]=C=C[CH][C]=C[O](22638)'],
    transitionState = 'TS61',
    kinetics = Arrhenius(A=(7.80481e+09,'s^-1'), n=1.09833, Ea=(177.686,'kJ/mol'), T0=(1,'K'), comment="""Estimated using an average for rate rule [R2H_S;Cd_rad_out_double;XH_out]
Euclidian distance = 0
Multiplied by reaction path degeneracy 2.0
family: intra_H_migration"""),
)

reaction(
    label = 'reaction62',
    reactants = ['[CH]=C=C[CH][C]=C[O](22638)'],
    products = ['[CH]=[C]C1C=[C]C1[O](26202)'],
    transitionState = 'TS62',
    kinetics = Arrhenius(A=(3.01156e+10,'s^-1'), n=0.428741, Ea=(355.449,'kJ/mol'), T0=(1,'K'), comment="""Estimated using average of templates [Rn;doublebond_intra_CdCdd;radadd_intra_cs] + [R4;doublebond_intra;radadd_intra_cs] for rate rule [R4;doublebond_intra_CdCdd;radadd_intra_cs]
Euclidian distance = 1.0
family: Intra_R_Add_Endocyclic"""),
)

reaction(
    label = 'reaction63',
    reactants = ['[CH]=C=C[CH][C]=C[O](22638)'],
    products = ['[O]C1[C]=CC=[C][CH]1(26203)'],
    transitionState = 'TS63',
    kinetics = Arrhenius(A=(2.16959e+10,'s^-1'), n=0.31, Ea=(35.7118,'kJ/mol'), T0=(1,'K'), comment="""Estimated using template [R6_linear;doublebond_intra;radadd_intra_cdsingleH] for rate rule [R6_linear;doublebond_intra_CdCdd;radadd_intra_cdsingleH]
Euclidian distance = 1.0
family: Intra_R_Add_Endocyclic
Ea raised from 32.1 to 35.7 kJ/mol to match endothermicity of reaction."""),
)

reaction(
    label = 'reaction64',
    reactants = ['[CH]=C=C[CH][C]=C[O](22638)'],
    products = ['[CH]=[C]CC=C=C=O(24786)'],
    transitionState = 'TS64',
    kinetics = Arrhenius(A=(5.14222e+08,'s^-1'), n=0.311, Ea=(34.991,'kJ/mol'), T0=(1,'K'), comment="""Estimated using template [R4;Y_rad;XH_Rrad] for rate rule [R4radExo;Y_rad;XH_Rrad]
Euclidian distance = 1.0
family: Intra_Disproportionation
Ea raised from 32.4 to 35.0 kJ/mol to match endothermicity of reaction."""),
)

reaction(
    label = 'reaction65',
    reactants = ['[CH]=C=C[C]=[C]C[O](26204)'],
    products = ['[CH]=C=C[CH][C]=C[O](22638)'],
    transitionState = 'TS65',
    kinetics = Arrhenius(A=(2.56742e+09,'s^-1'), n=1.0541, Ea=(193.078,'kJ/mol'), T0=(1,'K'), comment="""Estimated using template [R3Hall;Cd_rad_out_singleDe_Cd;XH_out] for rate rule [R3HJ;Cd_rad_out_singleDe_Cd;XH_out]
Euclidian distance = 1.0
Multiplied by reaction path degeneracy 2.0
family: intra_H_migration"""),
)

reaction(
    label = 'reaction66',
    reactants = ['[CH]=C=[C]C=[C]C[O](26205)'],
    products = ['[CH]=C=C[CH][C]=C[O](22638)'],
    transitionState = 'TS66',
    kinetics = Arrhenius(A=(2.26683e+09,'s^-1'), n=1.04717, Ea=(150.072,'kJ/mol'), T0=(1,'K'), comment="""Estimated using template [R4Hall;Cd_rad_out_double;XH_out] for rate rule [R4HJ_2;Cd_rad_out_double;XH_out]
Euclidian distance = 1.0
Multiplied by reaction path degeneracy 2.0
family: intra_H_migration"""),
)

reaction(
    label = 'reaction67',
    reactants = ['[CH]=[C]C[C][C]=C[O](26206)'],
    products = ['[CH]=C=C[CH][C]=C[O](22638)'],
    transitionState = 'TS67',
    kinetics = Arrhenius(A=(3.898e+11,'s^-1'), n=0.486, Ea=(22.8614,'kJ/mol'), T0=(1,'K'), comment="""Estimated using an average for rate rule [R2radExo;Y_rad;XH_Rrad]
Euclidian distance = 0
Multiplied by reaction path degeneracy 2.0
family: Intra_Disproportionation"""),
)

reaction(
    label = 'reaction68',
    reactants = ['[CH]=C=C[CH][C]=C[O](22638)'],
    products = ['C#CC=C[C]=C[O](22635)'],
    transitionState = 'TS68',
    kinetics = Arrhenius(A=(1e+10,'s^-1'), n=0, Ea=(0,'kJ/mol'), T0=(1,'K'), comment="""Estimated using an average for rate rule [NOS]
Euclidian distance = 0
family: 1,2-Birad_to_alkene"""),
)

reaction(
    label = 'reaction69',
    reactants = ['C#CC=O(21959)', '[CH][C]=[CH](21256)'],
    products = ['[CH]=C=C[CH][C]=C[O](22638)'],
    transitionState = 'TS69',
    kinetics = Arrhenius(A=(0.212494,'m^3/(mol*s)'), n=2.32278, Ea=(16.475,'kJ/mol'), T0=(1,'K'), comment="""Estimated using template [Ct-H_Ct-De;CJ] for rate rule [Ct-H_Ct-CO;CJ]
Euclidian distance = 1.0
Multiplied by reaction path degeneracy 2.0
family: R_Addition_MultipleBond"""),
)

reaction(
    label = 'reaction70',
    reactants = ['[CH]=O(373)', '[C][CH]C=C=[CH](21432)'],
    products = ['[CH]=C=C[CH][C]=C[O](22638)'],
    transitionState = 'TS70',
    kinetics = Arrhenius(A=(1.14854e+06,'m^3/(mol*s)'), n=0.575199, Ea=(34.3157,'kJ/mol'), T0=(1,'K'), comment="""Estimated using template [Y_rad;Birad] for rate rule [CO_pri_rad;Birad]
Euclidian distance = 2.0
family: Birad_R_Recombination"""),
)

reaction(
    label = 'reaction71',
    reactants = ['[CH]=C=C[CH][C]=C[O](22638)'],
    products = ['C#CC=[C]C=C[O](26207)'],
    transitionState = 'TS71',
    kinetics = Arrhenius(A=(1.949e+11,'s^-1'), n=0.486, Ea=(22.8614,'kJ/mol'), T0=(1,'K'), comment="""Estimated using an average for rate rule [R2radExo;Y_rad_De;XH_Rrad]
Euclidian distance = 0
family: Intra_Disproportionation"""),
)

reaction(
    label = 'reaction72',
    reactants = ['[CH]=C=C[CH][C]=C[O](22638)'],
    products = ['C#CC1[CH]C1=C[O](26162)'],
    transitionState = 'TS72',
    kinetics = Arrhenius(A=(2.18e+16,'s^-1'), n=0, Ea=(2.9288,'kJ/mol'), T0=(1,'K'), comment="""Estimated using template [R3;Y_rad_out;Ypri_rad_out] for rate rule [R3_SD;Y_rad_out;CdsinglepriDe_rad_out]
Euclidian distance = 2.2360679775
family: Birad_recombination"""),
)

reaction(
    label = 'reaction73',
    reactants = ['[CH]=C=C[CH][C]=C[O](22638)'],
    products = ['C#C[CH]C1[C]C1[O](26208)'],
    transitionState = 'TS73',
    kinetics = Arrhenius(A=(2.04797e+10,'s^-1'), n=0.473512, Ea=(252.024,'kJ/mol'), T0=(1,'K'), comment="""Estimated using template [R4;multiplebond_intra;radadd_intra_cs] for rate rule [R4;carbonylbond_intra_H;radadd_intra_cs]
Euclidian distance = 2.0
family: Intra_R_Add_Exocyclic"""),
)

reaction(
    label = 'reaction74',
    reactants = ['[C]#CC=C[C]C[O](26209)'],
    products = ['[CH]=C=C[CH][C]=C[O](22638)'],
    transitionState = 'TS74',
    kinetics = Arrhenius(A=(38203.4,'s^-1'), n=1.79759, Ea=(74.4436,'kJ/mol'), T0=(1,'K'), comment="""Estimated using template [R6Hall;Y_rad_out;XH_out] for rate rule [R6Hall;Ct_rad_out;XH_out]
Euclidian distance = 1.0
Multiplied by reaction path degeneracy 2.0
family: intra_H_migration"""),
)

reaction(
    label = 'reaction75',
    reactants = ['[C]#CC[CH][C][CH][O](26210)'],
    products = ['[CH]=C=C[CH][C]=C[O](22638)'],
    transitionState = 'TS75',
    kinetics = Arrhenius(A=(1.4874e+09,'s^-1'), n=1.045, Ea=(63.4002,'kJ/mol'), T0=(1,'K'), comment="""Estimated using an average for rate rule [R3radExo;Y_rad;XH_Rrad]
Euclidian distance = 0
Multiplied by reaction path degeneracy 2.0
family: Intra_Disproportionation"""),
)

reaction(
    label = 'reaction76',
    reactants = ['[C]#C[CH]C[C][CH][O](26211)'],
    products = ['[CH]=C=C[CH][C]=C[O](22638)'],
    transitionState = 'TS76',
    kinetics = Arrhenius(A=(1.02844e+09,'s^-1'), n=0.311, Ea=(24.9733,'kJ/mol'), T0=(1,'K'), comment="""Estimated using template [R4;Y_rad;XH_Rrad] for rate rule [R4radEndo;Y_rad;XH_Rrad]
Euclidian distance = 1.0
Multiplied by reaction path degeneracy 2.0
family: Intra_Disproportionation"""),
)

network(
    label = '4941',
    isomers = [
        '[CH]=C=C[CH][C]=C[O](22638)',
    ],
    reactants = [
        ('C#CC=O(21959)', '[CH]=C=[CH](18734)'),
    ],
    bathGas = {
        'N2': 0.25,
        'Ne': 0.25,
        'He': 0.25,
        'Ar': 0.25,
    },
)

pressureDependence(
    label = '4941',
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

