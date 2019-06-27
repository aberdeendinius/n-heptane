species(
    label = 'C#CC([O])C([O])C=C(22582)',
    structure = SMILES('C#CC([O])C([O])C=C'),
    E0 = (344.758,'kJ/mol'),
    modes = [
        HarmonicOscillator(frequencies=([2950,3100,1380,975,1025,1650,750,770,3400,2100,1380,1383.33,1386.67,1390,370,373.333,376.667,380,2800,3000,430,440,2175,525,3010,987.5,1337.5,450,1655,408.426,408.431,408.455,408.486],'cm^-1')),
        HinderedRotor(inertia=(0.191236,'amu*angstrom^2'), symmetry=1, barrier=(22.6485,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.191292,'amu*angstrom^2'), symmetry=1, barrier=(22.6488,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.191389,'amu*angstrom^2'), symmetry=1, barrier=(22.6493,'kJ/mol'), semiclassical=False),
    ],
    spinMultiplicity = 3,
    opticalIsomers = 1,
    molecularWeight = (110.111,'amu'),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[0.424544,0.0678809,-4.67261e-05,3.22755e-09,6.10177e-12,41603.2,31.6407], Tmin=(100,'K'), Tmax=(966.744,'K')), NASAPolynomial(coeffs=[18.6193,0.0169955,-5.62688e-06,9.90032e-10,-7.02363e-14,36945.2,-61.4194], Tmin=(966.744,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(344.758,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(320.107,'J/(mol*K)'), comment="""Thermo group additivity estimation: group(O2s-CsH) + group(O2s-CsH) + group(Cs-(Cds-Cds)CsOsH) + group(Cs-(Cds-Cds)CsOsH) + group(Cds-CdsCsH) + group(Cds-CdsHH) + group(Ct-CtCs) + group(Ct-CtH) + radical(CC(C)OJ) + radical(CC(C)OJ)"""),
)

species(
    label = 'C=CC=O(5269)',
    structure = SMILES('C=CC=O'),
    E0 = (-81.3387,'kJ/mol'),
    modes = [
        HarmonicOscillator(frequencies=([2782.5,750,1395,475,1775,1000,3010,987.5,1337.5,450,1655,2950,3100,1380,975,1025,1650],'cm^-1')),
        HinderedRotor(inertia=(0.873408,'amu*angstrom^2'), symmetry=1, barrier=(20.0814,'kJ/mol'), semiclassical=False),
    ],
    spinMultiplicity = 1,
    opticalIsomers = 1,
    molecularWeight = (56.0633,'amu'),
    collisionModel = TransportData(shapeIndex=2, epsilon=(3136.31,'J/mol'), sigma=(5.14154,'angstroms'), dipoleMoment=(0,'C*m'), polarizability=(0,'angstroms^3'), rotrelaxcollnum=0, comment="""Epsilon & sigma estimated with Tc=489.88 K, Pc=52.36 bar (from Joback method)"""),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[2.9738,0.0193269,-1.02836e-06,-7.40922e-09,2.6466e-12,-9743.32,12.1361], Tmin=(100,'K'), Tmax=(1315.19,'K')), NASAPolynomial(coeffs=[7.40832,0.0154746,-7.62321e-06,1.50372e-09,-1.06406e-13,-11743,-13.6408], Tmin=(1315.19,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(-81.3387,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(178.761,'J/(mol*K)'), comment="""Thermo group additivity estimation: group(Cd-Cd(CO)H) + group(Cds-O2d(Cds-Cds)H) + group(Cds-CdsHH)"""),
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
    label = 'C#CC([O])C1OC1[CH2](25919)',
    structure = SMILES('C#CC([O])C1OC1[CH2]'),
    E0 = (375.461,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
    molecularWeight = (110.111,'amu'),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[-0.729703,0.0823856,-8.69043e-05,4.57479e-08,-8.8636e-12,45346.5,31.2214], Tmin=(100,'K'), Tmax=(1522.96,'K')), NASAPolynomial(coeffs=[18.8867,0.0122625,4.8246e-07,-5.24907e-10,4.8741e-14,41528.7,-64.5857], Tmin=(1522.96,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(375.461,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(320.107,'J/(mol*K)'), comment="""Thermo group additivity estimation: group(O2s-CsCs) + group(O2s-CsH) + group(Cs-CsCsOsH) + group(Cs-CsCsOsH) + group(Cs-(Cds-Cds)CsOsH) + group(Cs-CsHHH) + group(Ct-CtCs) + group(Ct-CtH) + ring(Ethylene_oxide) + radical(CJCO) + radical(CC(C)OJ)"""),
)

species(
    label = '[CH]=C1OC1C([O])C=C(25920)',
    structure = SMILES('[CH]=C1OC1C([O])C=C'),
    E0 = (366.841,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
    molecularWeight = (110.111,'amu'),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[0.456343,0.0554554,1.03921e-05,-7.15411e-08,3.62158e-11,44269,28.3644], Tmin=(100,'K'), Tmax=(938.575,'K')), NASAPolynomial(coeffs=[25.7464,0.00456305,8.09713e-07,-1.56864e-10,8.16999e-16,37016,-105.391], Tmin=(938.575,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(366.841,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(324.264,'J/(mol*K)'), comment="""Thermo group additivity estimation: group(O2s-Cs(Cds-Cd)) + group(O2s-CsH) + group(Cs-(Cds-Cds)CsOsH) + group(Cs-(Cds-Cds)CsOsH) + group(Cds-CdsCsOs) + group(Cds-CdsCsH) + group(Cds-CdsHH) + group(Cds-CdsHH) + ring(methyleneoxirane) + radical(Cds_P) + radical(CC(C)OJ)"""),
)

species(
    label = 'C#CC1OC([CH2])C1[O](25450)',
    structure = SMILES('C#CC1OC([CH2])C1[O]'),
    E0 = (368.85,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
    molecularWeight = (110.111,'amu'),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[-0.665547,0.078499,-7.87499e-05,3.99139e-08,-7.45166e-12,44551.2,29.7934], Tmin=(100,'K'), Tmax=(1592.07,'K')), NASAPolynomial(coeffs=[17.6858,0.013683,-5.51515e-08,-4.20024e-10,4.09509e-14,41079,-59.8289], Tmin=(1592.07,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(368.85,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(324.264,'J/(mol*K)'), comment="""Thermo group additivity estimation: group(O2s-CsCs) + group(O2s-CsH) + group(Cs-CsCsOsH) + group(Cs-CsCsOsH) + group(Cs-(Cds-Cds)CsOsH) + group(Cs-CsHHH) + group(Ct-CtCs) + group(Ct-CtH) + ring(Oxetane) + radical(CC(C)OJ) + radical(CJC(C)OC)"""),
)

species(
    label = '[CH]=C1OC(C=C)C1[O](25491)',
    structure = SMILES('[CH]=C1OC(C=C)C1[O]'),
    E0 = (316.655,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
    molecularWeight = (110.111,'amu'),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[0.770177,0.0458528,4.08576e-05,-1.07047e-07,5.07585e-11,38224.6,25.5148], Tmin=(100,'K'), Tmax=(913.821,'K')), NASAPolynomial(coeffs=[25.8238,0.00229601,3.84e-06,-8.76007e-10,5.48345e-14,30885.5,-108.199], Tmin=(913.821,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(316.655,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(328.422,'J/(mol*K)'), comment="""Thermo group additivity estimation: group(O2s-Cs(Cds-Cd)) + group(O2s-CsH) + group(Cs-(Cds-Cds)CsOsH) + group(Cs-(Cds-Cds)CsOsH) + group(Cds-CdsCsOs) + group(Cds-CdsCsH) + group(Cds-CdsHH) + group(Cds-CdsHH) + ring(2methyleneoxetane) + radical(CC(C)OJ) + radical(Cds_P)"""),
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
    label = 'C#CC([O])C(=O)C=C(25921)',
    structure = SMILES('C#CC([O])C(=O)C=C'),
    E0 = (190.758,'kJ/mol'),
    modes = [
        HarmonicOscillator(frequencies=([750,770,3400,2100,2175,525,2950,3100,1380,975,1025,1650,1380,1390,370,380,2900,435,3010,987.5,1337.5,450,1655,375,552.5,462.5,1710,347.119,347.124,347.135],'cm^-1')),
        HinderedRotor(inertia=(0.00139911,'amu*angstrom^2'), symmetry=1, barrier=(0.119627,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.0620223,'amu*angstrom^2'), symmetry=1, barrier=(5.30241,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.339108,'amu*angstrom^2'), symmetry=1, barrier=(29.0031,'kJ/mol'), semiclassical=False),
    ],
    spinMultiplicity = 2,
    opticalIsomers = 1,
    molecularWeight = (109.103,'amu'),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[1.22884,0.0639264,-6.62453e-05,3.8447e-08,-9.30327e-12,23040.2,27.9585], Tmin=(100,'K'), Tmax=(983.819,'K')), NASAPolynomial(coeffs=[9.77721,0.0291706,-1.3254e-05,2.53843e-09,-1.78489e-13,21358.2,-13.143], Tmin=(983.819,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(190.758,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(295.164,'J/(mol*K)'), comment="""Thermo group additivity estimation: group(O2s-CsH) + group(Cs-(Cds-O2d)CsOsH) + group(Cds-O2d(Cds-Cds)Cs) + group(Cd-Cd(CO)H) + group(Cds-CdsHH) + group(Ct-CtCs) + group(Ct-CtH) + radical(C=OCOJ)"""),
)

species(
    label = 'C#CC(=O)C([O])C=C(25922)',
    structure = SMILES('C#CC(=O)C([O])C=C'),
    E0 = (190.441,'kJ/mol'),
    modes = [
        HarmonicOscillator(frequencies=([750,770,3400,2100,3010,987.5,1337.5,450,1655,2950,3100,1380,975,1025,1650,1380,1390,370,380,2900,435,2175,525,375,552.5,462.5,1710,180,180,4000],'cm^-1')),
        HinderedRotor(inertia=(0.224256,'amu*angstrom^2'), symmetry=1, barrier=(18.9682,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.823996,'amu*angstrom^2'), symmetry=1, barrier=(18.9453,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.0779199,'amu*angstrom^2'), symmetry=1, barrier=(45.7529,'kJ/mol'), semiclassical=False),
    ],
    spinMultiplicity = 2,
    opticalIsomers = 1,
    molecularWeight = (109.103,'amu'),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[1.23263,0.0629984,-6.29159e-05,3.45658e-08,-7.87749e-12,23002.7,28.0774], Tmin=(100,'K'), Tmax=(1042.95,'K')), NASAPolynomial(coeffs=[10.3534,0.0280182,-1.26069e-05,2.40803e-09,-1.69198e-13,21100.1,-16.3087], Tmin=(1042.95,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(190.441,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(295.164,'J/(mol*K)'), comment="""Thermo group additivity estimation: group(O2s-CsH) + group(Cs-(Cds-O2d)CsOsH) + group(Cds-CdsCsH) + group(Cds-O2d(Cds-Cds)Cs) + group(Cds-CdsHH) + group(Ct-CtCs) + group(Ct-CtH) + radical(C=OCOJ)"""),
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
    label = '[CH]=C(64)',
    structure = SMILES('[CH]=C'),
    E0 = (289.245,'kJ/mol'),
    modes = [
        HarmonicOscillator(frequencies=([2950,3100,1380,975,1025,1650,826.012,826.012,3240.27],'cm^-1')),
    ],
    spinMultiplicity = 2,
    opticalIsomers = 1,
    molecularWeight = (27.0452,'amu'),
    collisionModel = TransportData(shapeIndex=2, epsilon=(1737.73,'J/mol'), sigma=(4.1,'angstroms'), dipoleMoment=(0,'C*m'), polarizability=(0,'angstroms^3'), rotrelaxcollnum=1.0, comment="""GRI-Mech"""),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[3.90671,-0.00406241,3.8678e-05,-4.62976e-08,1.729e-11,34797.2,6.09789], Tmin=(100,'K'), Tmax=(931.962,'K')), NASAPolynomial(coeffs=[5.44797,0.00498356,-1.08821e-06,1.79837e-10,-1.45096e-14,33829.8,-4.87808], Tmin=(931.962,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(289.245,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(108.088,'J/(mol*K)'), comment="""Thermo group additivity estimation: group(Cds-CdsHH) + group(Cds-CdsHH) + radical(Cds_P)"""),
)

species(
    label = 'C#CC([O])C=O(22474)',
    structure = SMILES('C#CC([O])C=O'),
    E0 = (139.346,'kJ/mol'),
    modes = [
        HarmonicOscillator(frequencies=([2782.5,750,1395,475,1775,1000,750,770,3400,2100,1380,1390,370,380,2900,435,2175,525,474.687],'cm^-1')),
        HinderedRotor(inertia=(0.059994,'amu*angstrom^2'), symmetry=1, barrier=(9.59305,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.659976,'amu*angstrom^2'), symmetry=1, barrier=(105.493,'kJ/mol'), semiclassical=False),
    ],
    spinMultiplicity = 2,
    opticalIsomers = 1,
    molecularWeight = (83.0654,'amu'),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[1.87851,0.0410031,-3.92309e-05,1.95951e-08,-3.83768e-12,16840.6,22.268], Tmin=(100,'K'), Tmax=(1295.51,'K')), NASAPolynomial(coeffs=[11.1826,0.0111839,-3.44062e-06,5.26944e-10,-3.24794e-14,14521.5,-24.6743], Tmin=(1295.51,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(139.346,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(199.547,'J/(mol*K)'), comment="""Thermo group additivity estimation: group(O2s-CsH) + group(Cs-(Cds-O2d)CsOsH) + group(Cds-OdCsH) + group(Ct-CtCs) + group(Ct-CtH) + radical(C=OCOJ)"""),
)

species(
    label = '[CH2]C=C[O](5266)',
    structure = SMILES('[CH2]C=C[O]'),
    E0 = (90.2929,'kJ/mol'),
    modes = [
        HarmonicOscillator(frequencies=([2995,3025,975,1000,1300,1375,400,500,1630,1680,3000,3100,440,815,1455,1000,180],'cm^-1')),
        HinderedRotor(inertia=(1.57685,'amu*angstrom^2'), symmetry=1, barrier=(36.2549,'kJ/mol'), semiclassical=False),
    ],
    spinMultiplicity = 3,
    opticalIsomers = 1,
    molecularWeight = (56.0633,'amu'),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[2.69019,0.0144913,4.15491e-05,-7.27602e-08,3.14101e-11,10920.2,13.4175], Tmin=(100,'K'), Tmax=(922.751,'K')), NASAPolynomial(coeffs=[14.044,0.00224417,1.35973e-06,-3.04875e-10,1.62832e-14,7250.86,-48.974], Tmin=(922.751,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(90.2929,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(178.761,'J/(mol*K)'), comment="""Thermo group additivity estimation: group(O2s-(Cds-Cd)H) + group(Cs-(Cds-Cds)HHH) + group(Cds-CdsCsH) + group(Cds-CdsOsH) + radical(Allyl_P) + radical(C=COJ)"""),
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
    label = 'C=CC([O])C=O(12780)',
    structure = SMILES('C=CC([O])C=O'),
    E0 = (-26.6044,'kJ/mol'),
    modes = [
        HarmonicOscillator(frequencies=([2782.5,750,1395,475,1775,1000,2950,3100,1380,975,1025,1650,1380,1390,370,380,2900,435,3010,987.5,1337.5,450,1655,180,4000],'cm^-1')),
        HinderedRotor(inertia=(0.796353,'amu*angstrom^2'), symmetry=1, barrier=(18.3097,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.0379599,'amu*angstrom^2'), symmetry=1, barrier=(18.311,'kJ/mol'), semiclassical=False),
    ],
    spinMultiplicity = 2,
    opticalIsomers = 1,
    molecularWeight = (85.0813,'amu'),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[2.03604,0.0365188,-1.44093e-05,-5.26523e-09,3.93449e-12,-3123.41,24.0883], Tmin=(100,'K'), Tmax=(1084.25,'K')), NASAPolynomial(coeffs=[10.1946,0.0190919,-7.83047e-06,1.46841e-09,-1.03414e-13,-5637.41,-19.3665], Tmin=(1084.25,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(-26.6044,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(249.434,'J/(mol*K)'), comment="""Thermo group additivity estimation: group(O2s-CsH) + group(Cs-(Cds-O2d)CsOsH) + group(Cds-CdsCsH) + group(Cds-OdCsH) + group(Cds-CdsHH) + radical(C=OCOJ)"""),
)

species(
    label = 'C#CC([O])C(O)=C[CH2](25923)',
    structure = SMILES('C#CC([O])C(O)=C[CH2]'),
    E0 = (208.142,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
    molecularWeight = (110.111,'amu'),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[-0.729066,0.0845515,-8.84781e-05,4.42945e-08,-8.29441e-12,25220.7,32.4512], Tmin=(100,'K'), Tmax=(1489.35,'K')), NASAPolynomial(coeffs=[23.1633,0.0090008,-9.23449e-07,-2.82758e-11,6.89052e-15,19366.3,-88.0946], Tmin=(1489.35,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(208.142,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(315.95,'J/(mol*K)'), comment="""Thermo group additivity estimation: group(O2s-CsH) + group(O2s-(Cds-Cd)H) + group(Cs-(Cds-Cds)(Cds-Cds)OsH) + group(Cs-(Cds-Cds)HHH) + group(Cds-CdsCsOs) + group(Cds-CdsCsH) + group(Ct-CtCs) + group(Ct-CtH) + radical(CC(C)OJ) + radical(Allyl_P)"""),
)

species(
    label = '[CH]=C=C(O)C([O])C=C(25924)',
    structure = SMILES('[CH]=C=C(O)C([O])C=C'),
    E0 = (221.774,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
    molecularWeight = (110.111,'amu'),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[-0.435294,0.0823007,-8.64124e-05,4.39833e-08,-8.4763e-12,26846,33.1763], Tmin=(100,'K'), Tmax=(1400.9,'K')), NASAPolynomial(coeffs=[21.9931,0.011003,-2.29999e-06,2.57289e-10,-1.3142e-14,21274.2,-80.0469], Tmin=(1400.9,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(221.774,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(320.107,'J/(mol*K)'), comment="""Thermo group additivity estimation: group(O2s-CsH) + group(O2s-(Cds-Cd)H) + group(Cs-(Cds-Cds)(Cds-Cds)OsH) + group(Cds-CdsCsOs) + group(Cds-CdsCsH) + group(Cds-CdsHH) + group(Cds-CdsHH) + group(Cdd-CdsCds) + radical(C=C=CJ) + radical(CC(C)OJ)"""),
)

species(
    label = '[CH]=C=C([O])C(O)C=C(25925)',
    structure = SMILES('[CH]=C=C([O])C(O)C=C'),
    E0 = (129.218,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
    molecularWeight = (110.111,'amu'),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[0.244003,0.0768094,-8.10334e-05,4.24015e-08,-8.29039e-12,15681.7,32.4253], Tmin=(100,'K'), Tmax=(1001.95,'K')), NASAPolynomial(coeffs=[16.9918,0.0183934,-6.22273e-06,1.03698e-09,-6.83359e-14,11901.7,-50.5213], Tmin=(1001.95,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(129.218,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(320.107,'J/(mol*K)'), comment="""Thermo group additivity estimation: group(O2s-CsH) + group(O2s-(Cds-Cd)H) + group(Cs-(Cds-Cds)(Cds-Cds)OsH) + group(Cds-CdsCsOs) + group(Cds-CdsCsH) + group(Cds-CdsHH) + group(Cds-CdsHH) + group(Cdd-CdsCds) + radical(C=C(C)OJ) + radical(C=C=CJ)"""),
)

species(
    label = 'C#CC([O])C(O)[C]=C(25926)',
    structure = SMILES('C#CC([O])C(O)[C]=C'),
    E0 = (352.239,'kJ/mol'),
    modes = [
        HarmonicOscillator(frequencies=([750,770,3400,2100,2950,3100,1380,975,1025,1650,1380,1383.33,1386.67,1390,370,373.333,376.667,380,2800,3000,430,440,1685,370,2175,525,3615,1277.5,1000,180,398.586,632.254],'cm^-1')),
        HinderedRotor(inertia=(0.0635817,'amu*angstrom^2'), symmetry=1, barrier=(18.0346,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.0635804,'amu*angstrom^2'), symmetry=1, barrier=(18.0345,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.784373,'amu*angstrom^2'), symmetry=1, barrier=(18.0343,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.784388,'amu*angstrom^2'), symmetry=1, barrier=(18.0346,'kJ/mol'), semiclassical=False),
    ],
    spinMultiplicity = 3,
    opticalIsomers = 1,
    molecularWeight = (110.111,'amu'),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[0.101653,0.0779555,-8.24895e-05,4.38974e-08,-9.07601e-12,42511.6,33.4658], Tmin=(100,'K'), Tmax=(1190.89,'K')), NASAPolynomial(coeffs=[18.2285,0.0170699,-5.79931e-06,9.65267e-10,-6.32969e-14,38194.2,-57.1521], Tmin=(1190.89,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(352.239,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(315.95,'J/(mol*K)'), comment="""Thermo group additivity estimation: group(O2s-CsH) + group(O2s-CsH) + group(Cs-(Cds-Cds)CsOsH) + group(Cs-(Cds-Cds)CsOsH) + group(Cds-CdsCsH) + group(Cds-CdsHH) + group(Ct-CtCs) + group(Ct-CtH) + radical(CC(C)OJ) + radical(Cds_S)"""),
)

species(
    label = 'C#CC(O)C([O])=C[CH2](25927)',
    structure = SMILES('C#CC(O)C([O])=C[CH2]'),
    E0 = (115.586,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
    molecularWeight = (110.111,'amu'),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[-0.0996936,0.0796841,-8.53925e-05,4.57611e-08,-9.41417e-12,14058.4,31.8763], Tmin=(100,'K'), Tmax=(1269.1,'K')), NASAPolynomial(coeffs=[19.3292,0.0145759,-3.86289e-06,5.29387e-10,-3.04936e-14,9438.72,-65.2588], Tmin=(1269.1,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(115.586,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(315.95,'J/(mol*K)'), comment="""Thermo group additivity estimation: group(O2s-CsH) + group(O2s-(Cds-Cd)H) + group(Cs-(Cds-Cds)(Cds-Cds)OsH) + group(Cs-(Cds-Cds)HHH) + group(Cds-CdsCsOs) + group(Cds-CdsCsH) + group(Ct-CtCs) + group(Ct-CtH) + radical(Allyl_P) + radical(C=C(C)OJ)"""),
)

species(
    label = 'C#CC(O)C([O])[C]=C(25928)',
    structure = SMILES('C#CC(O)C([O])[C]=C'),
    E0 = (352.239,'kJ/mol'),
    modes = [
        HarmonicOscillator(frequencies=([750,770,3400,2100,2950,3100,1380,975,1025,1650,1380,1383.33,1386.67,1390,370,373.333,376.667,380,2800,3000,430,440,1685,370,2175,525,3615,1277.5,1000,180,398.586,632.254],'cm^-1')),
        HinderedRotor(inertia=(0.0635817,'amu*angstrom^2'), symmetry=1, barrier=(18.0346,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.0635804,'amu*angstrom^2'), symmetry=1, barrier=(18.0345,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.784373,'amu*angstrom^2'), symmetry=1, barrier=(18.0343,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.784388,'amu*angstrom^2'), symmetry=1, barrier=(18.0346,'kJ/mol'), semiclassical=False),
    ],
    spinMultiplicity = 3,
    opticalIsomers = 1,
    molecularWeight = (110.111,'amu'),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[0.101653,0.0779555,-8.24895e-05,4.38974e-08,-9.07601e-12,42511.6,33.4658], Tmin=(100,'K'), Tmax=(1190.89,'K')), NASAPolynomial(coeffs=[18.2285,0.0170699,-5.79931e-06,9.65267e-10,-6.32969e-14,38194.2,-57.1521], Tmin=(1190.89,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(352.239,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(315.95,'J/(mol*K)'), comment="""Thermo group additivity estimation: group(O2s-CsH) + group(O2s-CsH) + group(Cs-(Cds-Cds)CsOsH) + group(Cs-(Cds-Cds)CsOsH) + group(Cds-CdsCsH) + group(Cds-CdsHH) + group(Ct-CtCs) + group(Ct-CtH) + radical(CC(C)OJ) + radical(Cds_S)"""),
)

species(
    label = '[CH]=CC(O)C([O])C#C(25929)',
    structure = SMILES('[CH]=CC(O)C([O])C#C'),
    E0 = (361.494,'kJ/mol'),
    modes = [
        HarmonicOscillator(frequencies=([3120,650,792.5,1650,750,770,3400,2100,2175,525,1380,1383.33,1386.67,1390,370,373.333,376.667,380,2800,3000,430,440,3010,987.5,1337.5,450,1655,3615,1277.5,1000,370.042,370.042],'cm^-1')),
        HinderedRotor(inertia=(0.197524,'amu*angstrom^2'), symmetry=1, barrier=(19.1933,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.197524,'amu*angstrom^2'), symmetry=1, barrier=(19.1933,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.197524,'amu*angstrom^2'), symmetry=1, barrier=(19.1933,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.197524,'amu*angstrom^2'), symmetry=1, barrier=(19.1933,'kJ/mol'), semiclassical=False),
    ],
    spinMultiplicity = 3,
    opticalIsomers = 1,
    molecularWeight = (110.111,'amu'),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[-0.149888,0.0793161,-8.26679e-05,4.24944e-08,-8.38075e-12,43637.2,34.301], Tmin=(100,'K'), Tmax=(1313.12,'K')), NASAPolynomial(coeffs=[20.332,0.0136266,-3.86227e-06,5.72397e-10,-3.52323e-14,38542.5,-69.0092], Tmin=(1313.12,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(361.494,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(315.95,'J/(mol*K)'), comment="""Thermo group additivity estimation: group(O2s-CsH) + group(O2s-CsH) + group(Cs-(Cds-Cds)CsOsH) + group(Cs-(Cds-Cds)CsOsH) + group(Cds-CdsCsH) + group(Cds-CdsHH) + group(Ct-CtCs) + group(Ct-CtH) + radical(Cds_P) + radical(CC(C)OJ)"""),
)

species(
    label = '[C]#CC(O)C([O])C=C(25930)',
    structure = SMILES('[C]#CC(O)C([O])C=C'),
    E0 = (451.541,'kJ/mol'),
    modes = [
        HarmonicOscillator(frequencies=([2950,3100,1380,975,1025,1650,3010,987.5,1337.5,450,1655,3615,1277.5,1000,2175,525,1380,1383.33,1386.67,1390,370,373.333,376.667,380,2800,3000,430,440,452.361,452.413,452.415,452.428],'cm^-1')),
        HinderedRotor(inertia=(0.0906273,'amu*angstrom^2'), symmetry=1, barrier=(13.1692,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.0906451,'amu*angstrom^2'), symmetry=1, barrier=(13.1704,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.09066,'amu*angstrom^2'), symmetry=1, barrier=(13.1701,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.707343,'amu*angstrom^2'), symmetry=1, barrier=(102.74,'kJ/mol'), semiclassical=False),
    ],
    spinMultiplicity = 3,
    opticalIsomers = 1,
    molecularWeight = (110.111,'amu'),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[0.313023,0.0750827,-7.40407e-05,3.36543e-08,-4.67095e-12,54446,32.9752], Tmin=(100,'K'), Tmax=(947.169,'K')), NASAPolynomial(coeffs=[16.9546,0.0184591,-5.99365e-06,9.80607e-10,-6.44802e-14,50681,-49.6408], Tmin=(947.169,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(451.541,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(315.95,'J/(mol*K)'), comment="""Thermo group additivity estimation: group(O2s-CsH) + group(O2s-CsH) + group(Cs-(Cds-Cds)CsOsH) + group(Cs-(Cds-Cds)CsOsH) + group(Cds-CdsCsH) + group(Cds-CdsHH) + group(Ct-CtCs) + group(Ct-CtH) + radical(CC(C)OJ) + radical(Acetyl)"""),
)

species(
    label = '[CH]=CC([O])C(O)C#C(25931)',
    structure = SMILES('[CH]=CC([O])C(O)C#C'),
    E0 = (361.494,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
    molecularWeight = (110.111,'amu'),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[-0.149887,0.0793161,-8.26679e-05,4.24944e-08,-8.38074e-12,43637.2,34.301], Tmin=(100,'K'), Tmax=(1313.12,'K')), NASAPolynomial(coeffs=[20.332,0.0136266,-3.86227e-06,5.72399e-10,-3.52325e-14,38542.5,-69.0091], Tmin=(1313.12,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(361.494,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(315.95,'J/(mol*K)'), comment="""Thermo group additivity estimation: group(O2s-CsH) + group(O2s-CsH) + group(Cs-(Cds-Cds)CsOsH) + group(Cs-(Cds-Cds)CsOsH) + group(Cds-CdsCsH) + group(Cds-CdsHH) + group(Ct-CtCs) + group(Ct-CtH) + radical(CC(C)OJ) + radical(Cds_P)"""),
)

species(
    label = '[C]#CC([O])C(O)C=C(25932)',
    structure = SMILES('[C]#CC([O])C(O)C=C'),
    E0 = (451.541,'kJ/mol'),
    modes = [
        HarmonicOscillator(frequencies=([2950,3100,1380,975,1025,1650,3010,987.5,1337.5,450,1655,3615,1277.5,1000,2175,525,1380,1383.33,1386.67,1390,370,373.333,376.667,380,2800,3000,430,440,452.361,452.413,452.415,452.428],'cm^-1')),
        HinderedRotor(inertia=(0.0906273,'amu*angstrom^2'), symmetry=1, barrier=(13.1692,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.0906451,'amu*angstrom^2'), symmetry=1, barrier=(13.1704,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.09066,'amu*angstrom^2'), symmetry=1, barrier=(13.1701,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.707343,'amu*angstrom^2'), symmetry=1, barrier=(102.74,'kJ/mol'), semiclassical=False),
    ],
    spinMultiplicity = 3,
    opticalIsomers = 1,
    molecularWeight = (110.111,'amu'),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[0.313023,0.0750827,-7.40407e-05,3.36543e-08,-4.67095e-12,54446,32.9752], Tmin=(100,'K'), Tmax=(947.169,'K')), NASAPolynomial(coeffs=[16.9546,0.0184591,-5.99365e-06,9.80607e-10,-6.44802e-14,50681,-49.6408], Tmin=(947.169,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(451.541,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(315.95,'J/(mol*K)'), comment="""Thermo group additivity estimation: group(O2s-CsH) + group(O2s-CsH) + group(Cs-(Cds-Cds)CsOsH) + group(Cs-(Cds-Cds)CsOsH) + group(Cds-CdsCsH) + group(Cds-CdsHH) + group(Ct-CtCs) + group(Ct-CtH) + radical(Acetyl) + radical(CC(C)OJ)"""),
)

species(
    label = 'C#CC([O])[CH][O](23227)',
    structure = SMILES('C#CC([O])[CH][O]'),
    E0 = (453.158,'kJ/mol'),
    modes = [
        HarmonicOscillator(frequencies=([3025,407.5,1350,352.5,750,770,3400,2100,1380,1390,370,380,2900,435,2175,525,270.432,271.639,2343.48],'cm^-1')),
        HinderedRotor(inertia=(0.6648,'amu*angstrom^2'), symmetry=1, barrier=(35.16,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(1.41037,'amu*angstrom^2'), symmetry=1, barrier=(72.9656,'kJ/mol'), semiclassical=False),
    ],
    spinMultiplicity = 4,
    opticalIsomers = 1,
    molecularWeight = (83.0654,'amu'),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[1.52187,0.0587168,-9.13094e-05,7.51377e-08,-2.39962e-11,54587.5,22.7788], Tmin=(100,'K'), Tmax=(882.675,'K')), NASAPolynomial(coeffs=[8.70407,0.017867,-7.78132e-06,1.39446e-09,-9.1791e-14,53643,-9.14279], Tmin=(882.675,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(453.158,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(199.547,'J/(mol*K)'), comment="""Thermo group additivity estimation: group(O2s-CsH) + group(O2s-CsH) + group(Cs-(Cds-Cds)CsOsH) + group(Cs-CsOsHH) + group(Ct-CtCs) + group(Ct-CtH) + radical(CC(C)OJ) + radical(CCsJOH) + radical(CCOJ)"""),
)

species(
    label = 'C#CC([O])C([O])=C[CH2](25933)',
    structure = SMILES('C#CC([O])C([O])=C[CH2]'),
    E0 = (345.947,'kJ/mol'),
    modes = [
        HarmonicOscillator(frequencies=([750,770,3400,2100,350,440,435,1725,2175,525,3000,3100,440,815,1455,1000,3010,987.5,1337.5,450,1655,1380,1390,370,380,2900,435,411.071,412.854,412.9],'cm^-1')),
        HinderedRotor(inertia=(0.176703,'amu*angstrom^2'), symmetry=1, barrier=(21.4027,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.176738,'amu*angstrom^2'), symmetry=1, barrier=(21.4094,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.746751,'amu*angstrom^2'), symmetry=1, barrier=(90.0978,'kJ/mol'), semiclassical=False),
    ],
    spinMultiplicity = 4,
    opticalIsomers = 1,
    molecularWeight = (109.103,'amu'),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[0.021116,0.0761645,-8.10507e-05,4.26283e-08,-8.56831e-12,41760.8,31.2522], Tmin=(100,'K'), Tmax=(1313.7,'K')), NASAPolynomial(coeffs=[19.3689,0.0126089,-3.17867e-06,4.19047e-10,-2.36387e-14,37078.2,-65.8431], Tmin=(1313.7,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(345.947,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(295.164,'J/(mol*K)'), comment="""Thermo group additivity estimation: group(O2s-CsH) + group(O2s-(Cds-Cd)H) + group(Cs-(Cds-Cds)(Cds-Cds)OsH) + group(Cs-(Cds-Cds)HHH) + group(Cds-CdsCsOs) + group(Cds-CdsCsH) + group(Ct-CtCs) + group(Ct-CtH) + radical(Allyl_P) + radical(CC(C)OJ) + radical(C=C(C)OJ)"""),
)

species(
    label = 'C=CC([O])[CH][O](13593)',
    structure = SMILES('C=CC([O])[CH][O]'),
    E0 = (287.208,'kJ/mol'),
    modes = [
        HarmonicOscillator(frequencies=([3025,407.5,1350,352.5,2950,3100,1380,975,1025,1650,1380,1390,370,380,2900,435,3010,987.5,1337.5,450,1655,323.858,323.861,323.862,1611.9],'cm^-1')),
        HinderedRotor(inertia=(0.00160725,'amu*angstrom^2'), symmetry=1, barrier=(0.119627,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.185798,'amu*angstrom^2'), symmetry=1, barrier=(13.8288,'kJ/mol'), semiclassical=False),
    ],
    spinMultiplicity = 4,
    opticalIsomers = 1,
    molecularWeight = (85.0813,'amu'),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[1.7852,0.0528407,-6.0998e-05,4.24067e-08,-1.25977e-11,34619.2,24.2296], Tmin=(100,'K'), Tmax=(802.732,'K')), NASAPolynomial(coeffs=[6.94132,0.0271475,-1.29866e-05,2.53281e-09,-1.79327e-13,33791.4,0.487328], Tmin=(802.732,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(287.208,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(249.434,'J/(mol*K)'), comment="""Thermo group additivity estimation: group(O2s-CsH) + group(O2s-CsH) + group(Cs-(Cds-Cds)CsOsH) + group(Cs-CsOsHH) + group(Cds-CdsCsH) + group(Cds-CdsHH) + radical(CCsJOH) + radical(CC(C)OJ) + radical(CCOJ)"""),
)

species(
    label = '[CH]=C=C([O])C([O])C=C(25934)',
    structure = SMILES('[CH]=C=C([O])C([O])C=C'),
    E0 = (359.579,'kJ/mol'),
    modes = [
        HarmonicOscillator(frequencies=([3120,650,792.5,1650,540,610,2055,350,440,435,1725,2950,3100,1380,975,1025,1650,1380,1390,370,380,2900,435,3010,987.5,1337.5,450,1655,221.034,223.039,223.369],'cm^-1')),
        HinderedRotor(inertia=(3.41346,'amu*angstrom^2'), symmetry=1, barrier=(119.627,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.372005,'amu*angstrom^2'), symmetry=1, barrier=(12.9618,'kJ/mol'), semiclassical=False),
    ],
    spinMultiplicity = 4,
    opticalIsomers = 1,
    molecularWeight = (109.103,'amu'),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[0.41365,0.0727326,-7.48276e-05,3.69678e-08,-6.51178e-12,43381.9,31.6247], Tmin=(100,'K'), Tmax=(993.988,'K')), NASAPolynomial(coeffs=[16.8966,0.0166524,-5.66713e-06,9.56694e-10,-6.39514e-14,39598.8,-50.3438], Tmin=(993.988,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(359.579,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(299.321,'J/(mol*K)'), comment="""Thermo group additivity estimation: group(O2s-CsH) + group(O2s-(Cds-Cd)H) + group(Cs-(Cds-Cds)(Cds-Cds)OsH) + group(Cds-CdsCsOs) + group(Cds-CdsCsH) + group(Cds-CdsHH) + group(Cds-CdsHH) + group(Cdd-CdsCds) + radical(CC(C)OJ) + radical(C=C=CJ) + radical(C=C(C)OJ)"""),
)

species(
    label = 'C#CC([O])C([O])[C]=C(25935)',
    structure = SMILES('C#CC([O])C([O])[C]=C'),
    E0 = (582.6,'kJ/mol'),
    modes = [
        HarmonicOscillator(frequencies=([2950,3100,1380,975,1025,1650,750,770,3400,2100,1380,1383.33,1386.67,1390,370,373.333,376.667,380,2800,3000,430,440,2175,525,1685,370,403.444,403.465,403.518,403.52],'cm^-1')),
        HinderedRotor(inertia=(0.167387,'amu*angstrom^2'), symmetry=1, barrier=(19.3366,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.167344,'amu*angstrom^2'), symmetry=1, barrier=(19.3369,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(1.03551,'amu*angstrom^2'), symmetry=1, barrier=(119.627,'kJ/mol'), semiclassical=False),
    ],
    spinMultiplicity = 4,
    opticalIsomers = 1,
    molecularWeight = (109.103,'amu'),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[0.35868,0.0728693,-7.28722e-05,3.42453e-08,-5.59965e-12,70208,32.3506], Tmin=(100,'K'), Tmax=(1018.01,'K')), NASAPolynomial(coeffs=[17.6126,0.0161875,-5.728e-06,9.97527e-10,-6.81328e-14,66119.3,-54.0256], Tmin=(1018.01,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(582.6,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(295.164,'J/(mol*K)'), comment="""Thermo group additivity estimation: group(O2s-CsH) + group(O2s-CsH) + group(Cs-(Cds-Cds)CsOsH) + group(Cs-(Cds-Cds)CsOsH) + group(Cds-CdsCsH) + group(Cds-CdsHH) + group(Ct-CtCs) + group(Ct-CtH) + radical(CC(C)OJ) + radical(Cds_S) + radical(CC(C)OJ)"""),
)

species(
    label = '[CH]=CC([O])C([O])C#C(25936)',
    structure = SMILES('[CH]=CC([O])C([O])C#C'),
    E0 = (591.854,'kJ/mol'),
    modes = [
        HarmonicOscillator(frequencies=([3120,650,792.5,1650,750,770,3400,2100,1380,1383.33,1386.67,1390,370,373.333,376.667,380,2800,3000,430,440,2175,525,3010,987.5,1337.5,450,1655,357.875,357.932,358.52],'cm^-1')),
        HinderedRotor(inertia=(0.259364,'amu*angstrom^2'), symmetry=1, barrier=(23.6396,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.259866,'amu*angstrom^2'), symmetry=1, barrier=(23.6358,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.257529,'amu*angstrom^2'), symmetry=1, barrier=(23.6381,'kJ/mol'), semiclassical=False),
    ],
    spinMultiplicity = 4,
    opticalIsomers = 1,
    molecularWeight = (109.103,'amu'),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[-0.0376045,0.0758863,-7.85945e-05,3.96526e-08,-7.63738e-12,71340,33.7082], Tmin=(100,'K'), Tmax=(1364.52,'K')), NASAPolynomial(coeffs=[20.4074,0.0116036,-3.14767e-06,4.55221e-10,-2.78316e-14,66165.4,-69.7975], Tmin=(1364.52,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(591.854,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(295.164,'J/(mol*K)'), comment="""Thermo group additivity estimation: group(O2s-CsH) + group(O2s-CsH) + group(Cs-(Cds-Cds)CsOsH) + group(Cs-(Cds-Cds)CsOsH) + group(Cds-CdsCsH) + group(Cds-CdsHH) + group(Ct-CtCs) + group(Ct-CtH) + radical(CC(C)OJ) + radical(CC(C)OJ) + radical(Cds_P)"""),
)

species(
    label = '[C]#CC([O])C([O])C=C(25937)',
    structure = SMILES('[C]#CC([O])C([O])C=C'),
    E0 = (681.902,'kJ/mol'),
    modes = [
        HarmonicOscillator(frequencies=([2950,3100,1380,975,1025,1650,1380,1383.33,1386.67,1390,370,373.333,376.667,380,2800,3000,430,440,2175,525,3010,987.5,1337.5,450,1655,465.704,465.704,465.704,465.704,465.704],'cm^-1')),
        HinderedRotor(inertia=(0.10207,'amu*angstrom^2'), symmetry=1, barrier=(15.7088,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.583581,'amu*angstrom^2'), symmetry=1, barrier=(89.8149,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.10207,'amu*angstrom^2'), symmetry=1, barrier=(15.7088,'kJ/mol'), semiclassical=False),
    ],
    spinMultiplicity = 4,
    opticalIsomers = 1,
    molecularWeight = (109.103,'amu'),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[0.469905,0.0711606,-6.83941e-05,2.8969e-08,-3.22099e-12,82146.8,32.2201], Tmin=(100,'K'), Tmax=(952.003,'K')), NASAPolynomial(coeffs=[16.906,0.0166384,-5.39186e-06,8.89374e-10,-5.91858e-14,78358.6,-49.7256], Tmin=(952.003,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(681.902,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(295.164,'J/(mol*K)'), comment="""Thermo group additivity estimation: group(O2s-CsH) + group(O2s-CsH) + group(Cs-(Cds-Cds)CsOsH) + group(Cs-(Cds-Cds)CsOsH) + group(Cds-CdsCsH) + group(Cds-CdsHH) + group(Ct-CtCs) + group(Ct-CtH) + radical(CC(C)OJ) + radical(Acetyl) + radical(CC(C)OJ)"""),
)

species(
    label = 'C#CC([O])C1[CH]CO1(25938)',
    structure = SMILES('C#CC([O])C1[CH]CO1'),
    E0 = (370.885,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
    molecularWeight = (110.111,'amu'),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[1.13573,0.0478122,1.53284e-05,-6.99753e-08,3.67606e-11,44724.9,27.8109], Tmin=(100,'K'), Tmax=(869.573,'K')), NASAPolynomial(coeffs=[18.4545,0.0122048,7.50855e-07,-5.33368e-10,4.49272e-14,40047.2,-62.8998], Tmin=(869.573,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(370.885,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(324.264,'J/(mol*K)'), comment="""Thermo group additivity estimation: group(O2s-CsCs) + group(O2s-CsH) + group(Cs-CsCsOsH) + group(Cs-CsCsHH) + group(Cs-(Cds-Cds)CsOsH) + group(Cs-CsOsHH) + group(Ct-CtCs) + group(Ct-CtH) + ring(Oxetane) + radical(CC(C)OJ) + radical(CCJCO)"""),
)

species(
    label = 'C=CC([O])C1[C]=CO1(25859)',
    structure = SMILES('C=CC([O])C1[C]=CO1'),
    E0 = (338.121,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
    molecularWeight = (110.111,'amu'),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[0.632794,0.047374,3.95339e-05,-1.04428e-07,4.83955e-11,40812.3,28.3745], Tmin=(100,'K'), Tmax=(937.078,'K')), NASAPolynomial(coeffs=[26.9029,0.00308634,1.81959e-06,-3.30788e-10,1.02228e-14,32910,-112.551], Tmin=(937.078,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(338.121,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(324.264,'J/(mol*K)'), comment="""Thermo group additivity estimation: group(O2s-Cs(Cds-Cd)) + group(O2s-CsH) + group(Cs-(Cds-Cds)CsOsH) + group(Cs-(Cds-Cds)CsOsH) + group(Cds-CdsCsH) + group(Cds-CdsCsH) + group(Cds-CdsOsH) + group(Cds-CdsHH) + ring(Cyclobutene) + radical(CC(C)OJ) + radical(Cds_S)"""),
)

species(
    label = 'C#CC1OC[CH]C1[O](25774)',
    structure = SMILES('C#CC1OC[CH]C1[O]'),
    E0 = (292.423,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
    molecularWeight = (110.111,'amu'),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[1.47115,0.0362031,4.97138e-05,-1.05566e-07,4.91941e-11,35279.9,23.7121], Tmin=(100,'K'), Tmax=(882.198,'K')), NASAPolynomial(coeffs=[18.7489,0.0111891,1.57553e-06,-6.70623e-10,5.16843e-14,30156.3,-69.2387], Tmin=(882.198,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(292.423,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(328.422,'J/(mol*K)'), comment="""Thermo group additivity estimation: group(O2s-CsCs) + group(O2s-CsH) + group(Cs-CsCsOsH) + group(Cs-CsCsHH) + group(Cs-(Cds-Cds)CsOsH) + group(Cs-CsOsHH) + group(Ct-CtCs) + group(Ct-CtH) + ring(Tetrahydrofuran) + radical(CCJCO) + radical(CC(C)OJ)"""),
)

species(
    label = 'C=CC1OC=[C]C1[O](25541)',
    structure = SMILES('C=CC1OC=[C]C1[O]'),
    E0 = (234.115,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
    molecularWeight = (110.111,'amu'),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[1.19383,0.0353762,6.44965e-05,-1.27749e-07,5.73201e-11,28283.1,24.1019], Tmin=(100,'K'), Tmax=(913.931,'K')), NASAPolynomial(coeffs=[24.1644,0.00369316,3.49216e-06,-8.18421e-10,5.05031e-14,21208.9,-100.382], Tmin=(913.931,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(234.115,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(328.422,'J/(mol*K)'), comment="""Thermo group additivity estimation: group(O2s-Cs(Cds-Cd)) + group(O2s-CsH) + group(Cs-(Cds-Cds)CsOsH) + group(Cs-(Cds-Cds)CsOsH) + group(Cds-CdsCsH) + group(Cds-CdsCsH) + group(Cds-CdsOsH) + group(Cds-CdsHH) + ring(2,3-Dihydrofuran) + radical(Cds_S) + radical(CC(C)OJ)"""),
)

species(
    label = 'C#CC(O)C(=O)C=C(25939)',
    structure = SMILES('C#CC(O)C(=O)C=C'),
    E0 = (-52.9749,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
    molecularWeight = (110.111,'amu'),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[0.969854,0.0671557,-6.38741e-05,3.25524e-08,-6.8158e-12,-6262.71,28.1007], Tmin=(100,'K'), Tmax=(1134.58,'K')), NASAPolynomial(coeffs=[11.9876,0.0283116,-1.25186e-05,2.37605e-09,-1.66455e-13,-8762.79,-26.4446], Tmin=(1134.58,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(-52.9749,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(315.95,'J/(mol*K)'), comment="""Thermo group additivity estimation: group(O2s-CsH) + group(Cs-(Cds-O2d)CsOsH) + group(Cds-O2d(Cds-Cds)Cs) + group(Cd-Cd(CO)H) + group(Cds-CdsHH) + group(Ct-CtCs) + group(Ct-CtH)"""),
)

species(
    label = 'C#CC(=O)C(O)C=C(25940)',
    structure = SMILES('C#CC(=O)C(O)C=C'),
    E0 = (-53.2921,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
    molecularWeight = (110.111,'amu'),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[0.931571,0.0667038,-6.21098e-05,3.0558e-08,-6.13511e-12,-6298.45,28.3718], Tmin=(100,'K'), Tmax=(1184.79,'K')), NASAPolynomial(coeffs=[12.716,0.0269186,-1.17403e-05,2.21603e-09,-1.54793e-13,-9090.88,-30.4795], Tmin=(1184.79,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(-53.2921,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(315.95,'J/(mol*K)'), comment="""Thermo group additivity estimation: group(O2s-CsH) + group(Cs-(Cds-O2d)CsOsH) + group(Cds-CdsCsH) + group(Cds-O2d(Cds-Cds)Cs) + group(Cds-CdsHH) + group(Ct-CtCs) + group(Ct-CtH)"""),
)

species(
    label = 'C#CC([O])[C]([O])C[CH2](25941)',
    structure = SMILES('C#CC([O])[C]([O])C[CH2]'),
    E0 = (599.186,'kJ/mol'),
    modes = [
        HarmonicOscillator(frequencies=([360,370,350,750,770,3400,2100,2750,2850,1437.5,1250,1305,750,350,3000,3100,440,815,1455,1000,2175,525,1380,1390,370,380,2900,435,378.196,378.197,378.197,1046.95],'cm^-1')),
        HinderedRotor(inertia=(0.00117859,'amu*angstrom^2'), symmetry=1, barrier=(0.119627,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.182579,'amu*angstrom^2'), symmetry=1, barrier=(18.5316,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.0928415,'amu*angstrom^2'), symmetry=1, barrier=(9.42333,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.699323,'amu*angstrom^2'), symmetry=1, barrier=(70.9804,'kJ/mol'), semiclassical=False),
    ],
    spinMultiplicity = 5,
    opticalIsomers = 1,
    molecularWeight = (110.111,'amu'),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[0.14639,0.085958,-0.00010745,7.04445e-08,-1.82553e-11,72203.2,32.4731], Tmin=(100,'K'), Tmax=(946.348,'K')), NASAPolynomial(coeffs=[14.6893,0.0244873,-1.00148e-05,1.80417e-09,-1.22044e-13,69450.7,-36.886], Tmin=(946.348,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(599.186,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(315.95,'J/(mol*K)'), comment="""Thermo group additivity estimation: group(O2s-CsH) + group(O2s-CsH) + group(Cs-CsCsOsH) + group(Cs-CsCsHH) + group(Cs-(Cds-Cds)CsOsH) + group(Cs-CsHHH) + group(Ct-CtCs) + group(Ct-CtH) + radical(RCCJ) + radical(CC(C)OJ) + radical(CC(C)OJ) + radical(C2CsJOH)"""),
)

species(
    label = '[CH]=C=C([O])C([O])C[CH2](25942)',
    structure = SMILES('[CH]=C=C([O])C([O])C[CH2]'),
    E0 = (440.8,'kJ/mol'),
    modes = [
        HarmonicOscillator(frequencies=([3120,650,792.5,1650,540,610,2055,350,440,435,1725,2750,2850,1437.5,1250,1305,750,350,3000,3100,440,815,1455,1000,1380,1390,370,380,2900,435,405.983,405.985,405.987],'cm^-1')),
        HinderedRotor(inertia=(0.00102276,'amu*angstrom^2'), symmetry=1, barrier=(0.119627,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.115449,'amu*angstrom^2'), symmetry=1, barrier=(13.5031,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.115448,'amu*angstrom^2'), symmetry=1, barrier=(13.5031,'kJ/mol'), semiclassical=False),
    ],
    spinMultiplicity = 5,
    opticalIsomers = 1,
    molecularWeight = (110.111,'amu'),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[0.153609,0.0788508,-8.35231e-05,4.33342e-08,-8.23783e-12,53159.6,34.0331], Tmin=(100,'K'), Tmax=(981.46,'K')), NASAPolynomial(coeffs=[17.4123,0.018372,-6.16102e-06,1.02147e-09,-6.72228e-14,49296.9,-51.3266], Tmin=(981.46,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(440.8,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(320.107,'J/(mol*K)'), comment="""Thermo group additivity estimation: group(O2s-CsH) + group(O2s-(Cds-Cd)H) + group(Cs-(Cds-Cds)CsOsH) + group(Cs-CsCsHH) + group(Cs-CsHHH) + group(Cds-CdsCsOs) + group(Cds-CdsHH) + group(Cdd-CdsCds) + radical(RCCJ) + radical(C=C(C)OJ) + radical(C=C=CJ) + radical(CC(C)OJ)"""),
)

species(
    label = '[CH]=C=C([O])C([O])[CH]C(25943)',
    structure = SMILES('[CH]=C=C([O])C([O])[CH]C'),
    E0 = (435.456,'kJ/mol'),
    modes = [
        HarmonicOscillator(frequencies=([3120,650,792.5,1650,2750,2800,2850,1350,1500,750,1050,1375,1000,540,610,2055,350,440,435,1725,1380,1390,370,380,2900,435,3025,407.5,1350,352.5,439.904,439.905,439.906],'cm^-1')),
        HinderedRotor(inertia=(0.000871125,'amu*angstrom^2'), symmetry=1, barrier=(0.119627,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.0946991,'amu*angstrom^2'), symmetry=1, barrier=(13.0044,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.0946993,'amu*angstrom^2'), symmetry=1, barrier=(13.0043,'kJ/mol'), semiclassical=False),
    ],
    spinMultiplicity = 5,
    opticalIsomers = 1,
    molecularWeight = (110.111,'amu'),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[-0.0683534,0.0773273,-7.96833e-05,4.1012e-08,-8.09237e-12,52530.1,35.5422], Tmin=(100,'K'), Tmax=(1339.73,'K')), NASAPolynomial(coeffs=[19.4112,0.0143178,-3.70591e-06,5.02577e-10,-2.89012e-14,47745.8,-62.5082], Tmin=(1339.73,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(435.456,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(320.107,'J/(mol*K)'), comment="""Thermo group additivity estimation: group(O2s-CsH) + group(O2s-(Cds-Cd)H) + group(Cs-(Cds-Cds)CsOsH) + group(Cs-CsCsHH) + group(Cs-CsHHH) + group(Cds-CdsCsOs) + group(Cds-CdsHH) + group(Cdd-CdsCds) + radical(C=C=CJ) + radical(CCJCO) + radical(CC(C)OJ) + radical(C=C(C)OJ)"""),
)

species(
    label = '[C]#CC([O])C([O])C[CH2](25944)',
    structure = SMILES('[C]#CC([O])C([O])C[CH2]'),
    E0 = (759.702,'kJ/mol'),
    modes = [
        HarmonicOscillator(frequencies=([2750,2850,1437.5,1250,1305,750,350,1380,1383.33,1386.67,1390,370,373.333,376.667,380,2800,3000,430,440,2175,525,3000,3100,440,815,1455,1000,200,800,1066.67,1333.33,1600],'cm^-1')),
        HinderedRotor(inertia=(0.156089,'amu*angstrom^2'), symmetry=1, barrier=(3.5888,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.156089,'amu*angstrom^2'), symmetry=1, barrier=(3.5888,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.156089,'amu*angstrom^2'), symmetry=1, barrier=(3.5888,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.156089,'amu*angstrom^2'), symmetry=1, barrier=(3.5888,'kJ/mol'), semiclassical=False),
    ],
    spinMultiplicity = 5,
    opticalIsomers = 1,
    molecularWeight = (110.111,'amu'),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[0.230154,0.0819232,-9.72698e-05,6.1009e-08,-1.5098e-11,91507.7,33.4862], Tmin=(100,'K'), Tmax=(991.936,'K')), NASAPolynomial(coeffs=[14.7532,0.0233598,-8.71195e-06,1.49143e-09,-9.79269e-14,88626.5,-36.4618], Tmin=(991.936,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(759.702,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(315.95,'J/(mol*K)'), comment="""Thermo group additivity estimation: group(O2s-CsH) + group(O2s-CsH) + group(Cs-CsCsOsH) + group(Cs-CsCsHH) + group(Cs-(Cds-Cds)CsOsH) + group(Cs-CsHHH) + group(Ct-CtCs) + group(Ct-CtH) + radical(RCCJ) + radical(CC(C)OJ) + radical(Acetyl) + radical(CC(C)OJ)"""),
)

species(
    label = '[C]#CC([O])C([O])[CH]C(25945)',
    structure = SMILES('[C]#CC([O])C([O])[CH]C'),
    E0 = (754.358,'kJ/mol'),
    modes = [
        HarmonicOscillator(frequencies=([2750,2800,2850,1350,1500,750,1050,1375,1000,1380,1383.33,1386.67,1390,370,373.333,376.667,380,2800,3000,430,440,2175,525,3025,407.5,1350,352.5,200,800,1066.67,1333.33,1600],'cm^-1')),
        HinderedRotor(inertia=(0.156089,'amu*angstrom^2'), symmetry=1, barrier=(3.5888,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.156089,'amu*angstrom^2'), symmetry=1, barrier=(3.5888,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.156089,'amu*angstrom^2'), symmetry=1, barrier=(3.5888,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.156089,'amu*angstrom^2'), symmetry=1, barrier=(3.5888,'kJ/mol'), semiclassical=False),
    ],
    spinMultiplicity = 5,
    opticalIsomers = 1,
    molecularWeight = (110.111,'amu'),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[0.231125,0.0778641,-8.49859e-05,4.83126e-08,-1.07658e-11,90868.3,34.1892], Tmin=(100,'K'), Tmax=(1103.3,'K')), NASAPolynomial(coeffs=[16.0766,0.0204156,-6.88018e-06,1.11676e-09,-7.1375e-14,87371.9,-43.8137], Tmin=(1103.3,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(754.358,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(315.95,'J/(mol*K)'), comment="""Thermo group additivity estimation: group(O2s-CsH) + group(O2s-CsH) + group(Cs-CsCsOsH) + group(Cs-CsCsHH) + group(Cs-(Cds-Cds)CsOsH) + group(Cs-CsHHH) + group(Ct-CtCs) + group(Ct-CtH) + radical(CC(C)OJ) + radical(CC(C)OJ) + radical(CCJCO) + radical(Acetyl)"""),
)

species(
    label = 'C#CC1OOC1C=C(22598)',
    structure = SMILES('C#CC1OOC1C=C'),
    E0 = (290.685,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
    molecularWeight = (110.111,'amu'),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[1.0167,0.0461175,2.95624e-05,-9.02524e-08,4.48873e-11,35087.3,25.417], Tmin=(100,'K'), Tmax=(884.326,'K')), NASAPolynomial(coeffs=[21.4339,0.00785252,2.72579e-06,-8.59502e-10,6.37691e-14,29361.4,-82.5318], Tmin=(884.326,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(290.685,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(324.264,'J/(mol*K)'), comment="""Thermo group additivity estimation: group(O2s-OsCs) + group(O2s-OsCs) + group(Cs-(Cds-Cds)CsOsH) + group(Cs-(Cds-Cds)CsOsH) + group(Cds-CdsCsH) + group(Cds-CdsHH) + group(Ct-CtCs) + group(Ct-CtH) + ring(12dioxetane)"""),
)

species(
    label = '[O]C=C=CCC=C[O](22588)',
    structure = SMILES('[O]C=C=CCC=C[O]'),
    E0 = (94.6684,'kJ/mol'),
    modes = [
        HarmonicOscillator(frequencies=([2750,2850,1437.5,1250,1305,750,350,2995,3005,3015,3025,975,983.333,991.667,1000,1300,1325,1350,1375,400,433.333,466.667,500,1630,1646.67,1663.33,1680,540,610,2055,180,180,180,180],'cm^-1')),
        HinderedRotor(inertia=(1.24294,'amu*angstrom^2'), symmetry=1, barrier=(28.5775,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(1.24272,'amu*angstrom^2'), symmetry=1, barrier=(28.5727,'kJ/mol'), semiclassical=False),
    ],
    spinMultiplicity = 3,
    opticalIsomers = 1,
    molecularWeight = (110.111,'amu'),
    collisionModel = TransportData(shapeIndex=2, epsilon=(4556.79,'J/mol'), sigma=(7.12408,'angstroms'), dipoleMoment=(0,'C*m'), polarizability=(0,'angstroms^3'), rotrelaxcollnum=0, comment="""Epsilon & sigma estimated with Tc=711.76 K, Pc=28.6 bar (from Joback method)"""),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[0.371633,0.0609434,-8.25709e-06,-4.93002e-08,2.78157e-11,11534,29.2454], Tmin=(100,'K'), Tmax=(935.092,'K')), NASAPolynomial(coeffs=[23.8145,0.00848651,-8.25121e-07,9.47677e-11,-1.27022e-14,5058.87,-93.4597], Tmin=(935.092,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(94.6684,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(324.264,'J/(mol*K)'), comment="""Thermo group additivity estimation: group(O2s-(Cds-Cd)H) + group(O2s-(Cds-Cd)H) + group(Cs-(Cds-Cds)(Cds-Cds)HH) + group(Cds-CdsCsH) + group(Cds-CdsCsH) + group(Cds-CdsOsH) + group(Cds-CdsOsH) + group(Cdd-CdsCds) + radical(C=COJ) + radical(C=COJ)"""),
)

species(
    label = '[O]C1[C]=CC[CH]C1[O](25946)',
    structure = SMILES('[O]C1[C]=CC[CH]C1[O]'),
    E0 = (531.919,'kJ/mol'),
    modes = [
        HarmonicOscillator(frequencies=([2750,2830,2910,2990,3070,3150,900,940,980,1020,1060,1100,300,800,800,800,800,800,800,800,800,800,800,800,1600,1600,1600,1600,1600,1600,1600,1600,1600,1600,1600,1600],'cm^-1')),
    ],
    spinMultiplicity = 5,
    opticalIsomers = 1,
    molecularWeight = (110.111,'amu'),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[1.43352,0.0373994,3.76196e-05,-7.85556e-08,3.31018e-11,64084.6,27.2466], Tmin=(100,'K'), Tmax=(973.383,'K')), NASAPolynomial(coeffs=[17.147,0.0195961,-7.01759e-06,1.37827e-09,-1.06083e-13,58810,-59.519], Tmin=(973.383,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(531.919,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(332.579,'J/(mol*K)'), comment="""Thermo group additivity estimation: group(O2s-CsH) + group(O2s-CsH) + group(Cs-CsCsOsH) + group(Cs-(Cds-Cds)CsOsH) + group(Cs-CsCsHH) + group(Cs-(Cds-Cds)CsHH) + group(Cds-CdsCsH) + group(Cds-CdsCsH) + ring(Cyclohexene) + radical(Cds_S) + radical(CCJCO) + radical(CC(C)OJ) + radical(CC(C)OJ)"""),
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
    label = 'C#CC([O])[CH]C=C(23684)',
    structure = SMILES('C#CC([O])[CH]C=C'),
    E0 = (404.121,'kJ/mol'),
    modes = [
        HarmonicOscillator(frequencies=([750,770,3400,2100,3025,407.5,1350,352.5,3010,987.5,1337.5,450,1655,2950,3100,1380,975,1025,1650,1380,1390,370,380,2900,435,2175,525,180,180,180],'cm^-1')),
        HinderedRotor(inertia=(1.73654,'amu*angstrom^2'), symmetry=1, barrier=(39.9264,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(1.738,'amu*angstrom^2'), symmetry=1, barrier=(39.96,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(1.73677,'amu*angstrom^2'), symmetry=1, barrier=(39.9318,'kJ/mol'), semiclassical=False),
    ],
    spinMultiplicity = 3,
    opticalIsomers = 1,
    molecularWeight = (94.1112,'amu'),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[0.745204,0.0611962,-3.72531e-05,-3.11168e-09,7.77261e-12,48731.1,24.614], Tmin=(100,'K'), Tmax=(961.514,'K')), NASAPolynomial(coeffs=[17.2269,0.016742,-5.51716e-06,9.64232e-10,-6.81434e-14,44447,-60.0495], Tmin=(961.514,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(404.121,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(295.164,'J/(mol*K)'), comment="""Thermo group additivity estimation: group(O2s-CsH) + group(Cs-(Cds-Cds)CsHH) + group(Cs-(Cds-Cds)CsOsH) + group(Cds-CdsCsH) + group(Cds-CdsHH) + group(Ct-CtCs) + group(Ct-CtH) + radical(CC(C)OJ) + radical(C=CCJCO)"""),
)

species(
    label = '[CH]=C=CC([O])C=C(20057)',
    structure = SMILES('[CH]=C=CC([O])C=C'),
    E0 = (436.333,'kJ/mol'),
    modes = [
        HarmonicOscillator(frequencies=([3120,650,792.5,1650,2950,3100,1380,975,1025,1650,540,610,2055,1380,1390,370,380,2900,435,2995,3025,975,1000,1300,1375,400,500,1630,1680,180,180],'cm^-1')),
        HinderedRotor(inertia=(0.970908,'amu*angstrom^2'), symmetry=1, barrier=(22.3231,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.9706,'amu*angstrom^2'), symmetry=1, barrier=(22.316,'kJ/mol'), semiclassical=False),
    ],
    spinMultiplicity = 3,
    opticalIsomers = 1,
    molecularWeight = (94.1112,'amu'),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[1.0173,0.0582997,-4.56264e-05,1.57995e-08,-1.35038e-12,52592.3,27.902], Tmin=(100,'K'), Tmax=(1075.48,'K')), NASAPolynomial(coeffs=[13.8762,0.0208721,-7.92724e-06,1.42018e-09,-9.74795e-14,49225.1,-37.8661], Tmin=(1075.48,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(436.333,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(299.321,'J/(mol*K)'), comment="""Thermo group additivity estimation: group(O2s-CsH) + group(Cs-(Cds-Cds)(Cds-Cds)OsH) + group(Cds-CdsCsH) + group(Cds-CdsCsH) + group(Cds-CdsHH) + group(Cds-CdsHH) + group(Cdd-CdsCds) + radical(C=C=CJ) + radical(CC(C)OJ)"""),
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
    E0 = (344.758,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS2',
    E0 = (452.914,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS3',
    E0 = (429.204,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS4',
    E0 = (466.801,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS5',
    E0 = (466.801,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS6',
    E0 = (434.78,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS7',
    E0 = (434.881,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS8',
    E0 = (344.758,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS9',
    E0 = (486.815,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS10',
    E0 = (344.758,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS11',
    E0 = (583.555,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS12',
    E0 = (507.33,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS13',
    E0 = (498.666,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS14',
    E0 = (420.039,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS15',
    E0 = (494.216,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS16',
    E0 = (420.039,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS17',
    E0 = (453.417,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS18',
    E0 = (405.802,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS19',
    E0 = (535.43,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS20',
    E0 = (456.203,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS21',
    E0 = (495.839,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS22',
    E0 = (360.168,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS23',
    E0 = (742.403,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS24',
    E0 = (557.752,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS25',
    E0 = (839.144,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS26',
    E0 = (571.384,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS27',
    E0 = (794.405,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS28',
    E0 = (803.659,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS29',
    E0 = (893.707,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS30',
    E0 = (477.269,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS31',
    E0 = (586.928,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS32',
    E0 = (400.795,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS33',
    E0 = (403.334,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS34',
    E0 = (408.158,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS35',
    E0 = (408.158,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS36',
    E0 = (621.487,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS37',
    E0 = (529.769,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS38',
    E0 = (474.681,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS39',
    E0 = (784.676,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS40',
    E0 = (779.331,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS41',
    E0 = (353.043,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS42',
    E0 = (545,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS43',
    E0 = (531.919,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS44',
    E0 = (811.002,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS45',
    E0 = (843.214,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

reaction(
    label = 'reaction1',
    reactants = ['C#CC([O])C([O])C=C(22582)'],
    products = ['C=CC=O(5269)', 'C#CC=O(21959)'],
    transitionState = 'TS1',
    kinetics = Arrhenius(A=(5e+12,'s^-1'), n=0, Ea=(0,'kJ/mol'), T0=(1,'K'), Tmin=(300,'K'), Tmax=(1500,'K'), comment="""Exact match found for rate rule [RJJ]
Euclidian distance = 0
family: 1,4_Linear_birad_scission"""),
)

reaction(
    label = 'reaction2',
    reactants = ['C#CC([O])C([O])C=C(22582)'],
    products = ['C#CC([O])C1OC1[CH2](25919)'],
    transitionState = 'TS2',
    kinetics = Arrhenius(A=(1e+10,'s^-1'), n=0, Ea=(108.156,'kJ/mol'), T0=(1,'K'), Tmin=(300,'K'), Tmax=(1500,'K'), comment="""Estimated using template [R4_S_D;doublebond_intra_2H_pri;radadd_intra] for rate rule [R4_S_D;doublebond_intra_2H_pri;radadd_intra_O]
Euclidian distance = 1.0
family: Intra_R_Add_Exocyclic"""),
)

reaction(
    label = 'reaction3',
    reactants = ['C#CC([O])C([O])C=C(22582)'],
    products = ['[CH]=C1OC1C([O])C=C(25920)'],
    transitionState = 'TS3',
    kinetics = Arrhenius(A=(3.98513e+09,'s^-1'), n=0.768, Ea=(84.4457,'kJ/mol'), T0=(1,'K'), comment="""Estimated using template [R4_S_T;triplebond_intra_H;radadd_intra] for rate rule [R4_S_T;triplebond_intra_H;radadd_intra_O]
Euclidian distance = 1.0
family: Intra_R_Add_Exocyclic"""),
)

reaction(
    label = 'reaction4',
    reactants = ['C#CC([O])C([O])C=C(22582)'],
    products = ['C#CC1OC([CH2])C1[O](25450)'],
    transitionState = 'TS4',
    kinetics = Arrhenius(A=(2.724e+10,'s^-1','*|/',3), n=0.478, Ea=(122.043,'kJ/mol'), T0=(1,'K'), Tmin=(600,'K'), Tmax=(2000,'K'), comment="""From training reaction 33 used for R5_SS_D;doublebond_intra_2H_pri;radadd_intra_O
Exact match found for rate rule [R5_SS_D;doublebond_intra_2H_pri;radadd_intra_O]
Euclidian distance = 0
family: Intra_R_Add_Exocyclic"""),
)

reaction(
    label = 'reaction5',
    reactants = ['C#CC([O])C([O])C=C(22582)'],
    products = ['[CH]=C1OC(C=C)C1[O](25491)'],
    transitionState = 'TS5',
    kinetics = Arrhenius(A=(2.724e+10,'s^-1'), n=0.478, Ea=(122.043,'kJ/mol'), T0=(1,'K'), comment="""Estimated using template [R5_SS;multiplebond_intra;radadd_intra_O] for rate rule [R5_SS_T;triplebond_intra_H;radadd_intra_O]
Euclidian distance = 2.2360679775
family: Intra_R_Add_Exocyclic"""),
)

reaction(
    label = 'reaction6',
    reactants = ['H(8)', 'C#CC([O])C(=O)C=C(25921)'],
    products = ['C#CC([O])C([O])C=C(22582)'],
    transitionState = 'TS6',
    kinetics = Arrhenius(A=(3.97e+07,'cm^3/(mol*s)'), n=1.88, Ea=(32.2168,'kJ/mol'), T0=(1,'K'), Tmin=(300,'K'), Tmax=(1500,'K'), comment="""From training reaction 2826 used for CO-CdCs_O;HJ
Exact match found for rate rule [CO-CdCs_O;HJ]
Euclidian distance = 0
family: R_Addition_MultipleBond"""),
)

reaction(
    label = 'reaction7',
    reactants = ['H(8)', 'C#CC(=O)C([O])C=C(25922)'],
    products = ['C#CC([O])C([O])C=C(22582)'],
    transitionState = 'TS7',
    kinetics = Arrhenius(A=(1.59e+07,'cm^3/(mol*s)'), n=1.84, Ea=(32.6352,'kJ/mol'), T0=(1,'K'), Tmin=(300,'K'), Tmax=(1500,'K'), comment="""From training reaction 2830 used for CO-CtCs_O;HJ
Exact match found for rate rule [CO-CtCs_O;HJ]
Euclidian distance = 0
family: R_Addition_MultipleBond"""),
)

reaction(
    label = 'reaction8',
    reactants = ['C=CC=O(5269)', '[CH]=C=C[O](8556)'],
    products = ['C#CC([O])C([O])C=C(22582)'],
    transitionState = 'TS8',
    kinetics = Arrhenius(A=(7.5e+06,'cm^3/(mol*s)'), n=2.16, Ea=(156.222,'kJ/mol'), T0=(1,'K'), Tmin=(300,'K'), Tmax=(1500,'K'), comment="""Estimated using template [CO-CdH_O;YJ] for rate rule [CO-CdH_O;CJ]
Euclidian distance = 1.0
family: R_Addition_MultipleBond
Ea raised from 151.2 to 156.2 kJ/mol to match endothermicity of reaction."""),
)

reaction(
    label = 'reaction9',
    reactants = ['[CH]=C(64)', 'C#CC([O])C=O(22474)'],
    products = ['C#CC([O])C([O])C=C(22582)'],
    transitionState = 'TS9',
    kinetics = Arrhenius(A=(2339.95,'m^3/(mol*s)'), n=0.573452, Ea=(58.2237,'kJ/mol'), T0=(1,'K'), comment="""Estimated using template [CO_O;CJ] for rate rule [CO-CsH_O;CdsJ-H]
Euclidian distance = 2.82842712475
family: R_Addition_MultipleBond"""),
)

reaction(
    label = 'reaction10',
    reactants = ['[CH2]C=C[O](5266)', 'C#CC=O(21959)'],
    products = ['C#CC([O])C([O])C=C(22582)'],
    transitionState = 'TS10',
    kinetics = Arrhenius(A=(2.99e+06,'cm^3/(mol*s)'), n=2.12, Ea=(170.171,'kJ/mol'), T0=(1,'K'), Tmin=(300,'K'), Tmax=(1500,'K'), comment="""Estimated using template [CO-CtH_O;YJ] for rate rule [CO-CtH_O;CJ]
Euclidian distance = 1.0
family: R_Addition_MultipleBond
Ea raised from 165.2 to 170.2 kJ/mol to match endothermicity of reaction."""),
)

reaction(
    label = 'reaction11',
    reactants = ['[C]#C(5143)', 'C=CC([O])C=O(12780)'],
    products = ['C#CC([O])C([O])C=C(22582)'],
    transitionState = 'TS11',
    kinetics = Arrhenius(A=(2339.95,'m^3/(mol*s)'), n=0.573452, Ea=(58.2237,'kJ/mol'), T0=(1,'K'), comment="""Estimated using template [CO_O;CJ] for rate rule [CO-CsH_O;CtJ_Ct]
Euclidian distance = 2.82842712475
family: R_Addition_MultipleBond"""),
)

reaction(
    label = 'reaction12',
    reactants = ['C#CC([O])C([O])C=C(22582)'],
    products = ['C#CC([O])C(O)=C[CH2](25923)'],
    transitionState = 'TS12',
    kinetics = Arrhenius(A=(1.52488e+09,'s^-1'), n=1.21745, Ea=(162.572,'kJ/mol'), T0=(1,'K'), comment="""Estimated using template [R2H_S;Y_rad_out;Cs_H_out_Cd] for rate rule [R2H_S;O_rad_out;Cs_H_out_Cd]
Euclidian distance = 1.0
family: intra_H_migration"""),
)

reaction(
    label = 'reaction13',
    reactants = ['C#CC([O])C([O])C=C(22582)'],
    products = ['[CH]=C=C(O)C([O])C=C(25924)'],
    transitionState = 'TS13',
    kinetics = Arrhenius(A=(1.70223e+09,'s^-1'), n=1.15155, Ea=(153.908,'kJ/mol'), T0=(1,'K'), comment="""Estimated using template [R2H_S;Y_rad_out;Cs_H_out_OneDe] for rate rule [R2H_S;O_rad_out;Cs_H_out_Ct]
Euclidian distance = 1.41421356237
family: intra_H_migration"""),
)

reaction(
    label = 'reaction14',
    reactants = ['C#CC([O])C([O])C=C(22582)'],
    products = ['[CH]=C=C([O])C(O)C=C(25925)'],
    transitionState = 'TS14',
    kinetics = Arrhenius(A=(111914,'s^-1'), n=2.27675, Ea=(75.2806,'kJ/mol'), T0=(1,'K'), comment="""Estimated using an average for rate rule [R3H_SS_Cs;O_rad_out;XH_out]
Euclidian distance = 0
family: intra_H_migration"""),
)

reaction(
    label = 'reaction15',
    reactants = ['C#CC([O])C(O)[C]=C(25926)'],
    products = ['C#CC([O])C([O])C=C(22582)'],
    transitionState = 'TS15',
    kinetics = Arrhenius(A=(2.4115e+09,'s^-1'), n=1.00333, Ea=(141.977,'kJ/mol'), T0=(1,'K'), comment="""Estimated using template [R3H_SS_Cs;Cd_rad_out_Cd;XH_out] for rate rule [R3H_SS_Cs;Cd_rad_out_Cd;O_H_out]
Euclidian distance = 1.0
family: intra_H_migration"""),
)

reaction(
    label = 'reaction16',
    reactants = ['C#CC([O])C([O])C=C(22582)'],
    products = ['C#CC(O)C([O])=C[CH2](25927)'],
    transitionState = 'TS16',
    kinetics = Arrhenius(A=(111914,'s^-1'), n=2.27675, Ea=(75.2806,'kJ/mol'), T0=(1,'K'), comment="""Estimated using an average for rate rule [R3H_SS_Cs;O_rad_out;XH_out]
Euclidian distance = 0
family: intra_H_migration"""),
)

reaction(
    label = 'reaction17',
    reactants = ['C#CC(O)C([O])[C]=C(25928)'],
    products = ['C#CC([O])C([O])C=C(22582)'],
    transitionState = 'TS17',
    kinetics = Arrhenius(A=(1.286e+08,'s^-1'), n=1.323, Ea=(101.177,'kJ/mol'), T0=(1,'K'), comment="""Estimated using template [R4H_RSR;Cd_rad_out_Cd;XH_out] for rate rule [R4H_SSS;Cd_rad_out_Cd;O_H_out]
Euclidian distance = 2.2360679775
family: intra_H_migration"""),
)

reaction(
    label = 'reaction18',
    reactants = ['[CH]=CC(O)C([O])C#C(25929)'],
    products = ['C#CC([O])C([O])C=C(22582)'],
    transitionState = 'TS18',
    kinetics = Arrhenius(A=(37100,'s^-1'), n=2.23, Ea=(44.3086,'kJ/mol'), T0=(1,'K'), comment="""Estimated using template [R4H_DSS;Cd_rad_out_singleH;XH_out] for rate rule [R4H_DSS;Cd_rad_out_singleH;O_H_out]
Euclidian distance = 1.0
family: intra_H_migration"""),
)

reaction(
    label = 'reaction19',
    reactants = ['[C]#CC(O)C([O])C=C(25930)'],
    products = ['C#CC([O])C([O])C=C(22582)'],
    transitionState = 'TS19',
    kinetics = Arrhenius(A=(0.00504539,'s^-1'), n=3.83, Ea=(83.8892,'kJ/mol'), T0=(1,'K'), comment="""Estimated using template [R4H_RSS;Y_rad_out;O_H_out] for rate rule [R4H_TSS;Ct_rad_out;O_H_out]
Euclidian distance = 1.41421356237
family: intra_H_migration"""),
)

reaction(
    label = 'reaction20',
    reactants = ['C#CC([O])C([O])C=C(22582)'],
    products = ['[CH]=CC([O])C(O)C#C(25931)'],
    transitionState = 'TS20',
    kinetics = Arrhenius(A=(2.468e+06,'s^-1','*|/',3), n=1.554, Ea=(111.445,'kJ/mol'), T0=(1,'K'), Tmin=(600,'K'), Tmax=(2000,'K'), comment="""From training reaction 290 used for R5H_SSSD;O_rad_out;Cd_H_out_singleH
Exact match found for rate rule [R5H_SSSD;O_rad_out;Cd_H_out_singleH]
Euclidian distance = 0
Multiplied by reaction path degeneracy 2.0
family: intra_H_migration"""),
)

reaction(
    label = 'reaction21',
    reactants = ['[C]#CC([O])C(O)C=C(25932)'],
    products = ['C#CC([O])C([O])C=C(22582)'],
    transitionState = 'TS21',
    kinetics = Arrhenius(A=(380071,'s^-1'), n=1.62386, Ea=(44.2978,'kJ/mol'), T0=(1,'K'), comment="""Estimated using template [R5H_RSSR;Y_rad_out;XH_out] for rate rule [R5H_TSSS;Ct_rad_out;O_H_out]
Euclidian distance = 2.44948974278
family: intra_H_migration"""),
)

reaction(
    label = 'reaction22',
    reactants = ['[CH2]C=C[O](5266)', '[CH]=C=C[O](8556)'],
    products = ['C#CC([O])C([O])C=C(22582)'],
    transitionState = 'TS22',
    kinetics = Arrhenius(A=(1.9789e+07,'m^3/(mol*s)'), n=-0.126319, Ea=(0,'kJ/mol'), T0=(1,'K'), comment="""Estimated using an average for rate rule [Y_rad;Y_rad]
Euclidian distance = 0
family: R_Recombination
Ea raised from -15.6 to -15.6 kJ/mol.
Ea raised from -15.6 to 0 kJ/mol."""),
)

reaction(
    label = 'reaction23',
    reactants = ['[CH]=C(64)', 'C#CC([O])[CH][O](23227)'],
    products = ['C#CC([O])C([O])C=C(22582)'],
    transitionState = 'TS23',
    kinetics = Arrhenius(A=(3.88428e+07,'m^3/(mol*s)'), n=0, Ea=(0,'kJ/mol'), T0=(1,'K'), comment="""Estimated using an average for rate rule [Y_rad;Cd_pri_rad]
Euclidian distance = 0
family: R_Recombination"""),
)

reaction(
    label = 'reaction24',
    reactants = ['H(8)', 'C#CC([O])C([O])=C[CH2](25933)'],
    products = ['C#CC([O])C([O])C=C(22582)'],
    transitionState = 'TS24',
    kinetics = Arrhenius(A=(4.34078e+06,'m^3/(mol*s)'), n=0.278577, Ea=(0,'kJ/mol'), T0=(1,'K'), comment="""Estimated using an average for rate rule [Y_rad;H_rad]
Euclidian distance = 0
family: R_Recombination
Ea raised from -1.4 to 0 kJ/mol."""),
)

reaction(
    label = 'reaction25',
    reactants = ['[C]#C(5143)', 'C=CC([O])[CH][O](13593)'],
    products = ['C#CC([O])C([O])C=C(22582)'],
    transitionState = 'TS25',
    kinetics = Arrhenius(A=(1.34536e+08,'m^3/(mol*s)'), n=0, Ea=(0,'kJ/mol'), T0=(1,'K'), comment="""Estimated using an average for rate rule [Ct_rad/Ct;Y_rad]
Euclidian distance = 0
family: R_Recombination"""),
)

reaction(
    label = 'reaction26',
    reactants = ['H(8)', '[CH]=C=C([O])C([O])C=C(25934)'],
    products = ['C#CC([O])C([O])C=C(22582)'],
    transitionState = 'TS26',
    kinetics = Arrhenius(A=(4.34078e+06,'m^3/(mol*s)'), n=0.278577, Ea=(0,'kJ/mol'), T0=(1,'K'), comment="""Estimated using an average for rate rule [Y_rad;H_rad]
Euclidian distance = 0
family: R_Recombination
Ea raised from -1.4 to 0 kJ/mol."""),
)

reaction(
    label = 'reaction27',
    reactants = ['H(8)', 'C#CC([O])C([O])[C]=C(25935)'],
    products = ['C#CC([O])C([O])C=C(22582)'],
    transitionState = 'TS27',
    kinetics = Arrhenius(A=(1e+13,'cm^3/(mol*s)'), n=0, Ea=(0,'kJ/mol'), T0=(1,'K'), comment="""From training reaction 40 used for Cd_rad/NonDe;H_rad
Exact match found for rate rule [Cd_rad/NonDe;H_rad]
Euclidian distance = 0
family: R_Recombination"""),
)

reaction(
    label = 'reaction28',
    reactants = ['H(8)', '[CH]=CC([O])C([O])C#C(25936)'],
    products = ['C#CC([O])C([O])C=C(22582)'],
    transitionState = 'TS28',
    kinetics = Arrhenius(A=(1.21e+14,'cm^3/(mol*s)','+|-',4.82e+13), n=0, Ea=(0,'kJ/mol'), T0=(1,'K'), Tmin=(298,'K'), comment="""From training reaction 60 used for H_rad;Cd_pri_rad
Exact match found for rate rule [Cd_pri_rad;H_rad]
Euclidian distance = 0
family: R_Recombination"""),
)

reaction(
    label = 'reaction29',
    reactants = ['H(8)', '[C]#CC([O])C([O])C=C(25937)'],
    products = ['C#CC([O])C([O])C=C(22582)'],
    transitionState = 'TS29',
    kinetics = Arrhenius(A=(1.81e+14,'cm^3/(mol*s)','*|/',3), n=0, Ea=(0,'kJ/mol'), T0=(1,'K'), Tmin=(300,'K'), Tmax=(2500,'K'), comment="""From training reaction 61 used for H_rad;Ct_rad/Ct
Exact match found for rate rule [H_rad;Ct_rad/Ct]
Euclidian distance = 0
family: R_Recombination"""),
)

reaction(
    label = 'reaction30',
    reactants = ['C#CC([O])C([O])C=C(22582)'],
    products = ['C#CC([O])C1[CH]CO1(25938)'],
    transitionState = 'TS30',
    kinetics = Arrhenius(A=(1.15968e+08,'s^-1'), n=1.10215, Ea=(132.51,'kJ/mol'), T0=(1,'K'), comment="""Estimated using template [R4_S_D;doublebond_intra_pri_2H;radadd_intra] for rate rule [R4_S_D;doublebond_intra_pri_2H;radadd_intra_O]
Euclidian distance = 1.0
family: Intra_R_Add_Endocyclic"""),
)

reaction(
    label = 'reaction31',
    reactants = ['C#CC([O])C([O])C=C(22582)'],
    products = ['C=CC([O])C1[C]=CO1(25859)'],
    transitionState = 'TS31',
    kinetics = Arrhenius(A=(1.503e+11,'s^-1'), n=0.221, Ea=(242.17,'kJ/mol'), T0=(1,'K'), comment="""Estimated using template [R4_S;multiplebond_intra;radadd_intra_O] for rate rule [R4_S_T;triplebond_intra_H;radadd_intra_O]
Euclidian distance = 2.2360679775
family: Intra_R_Add_Endocyclic"""),
)

reaction(
    label = 'reaction32',
    reactants = ['C#CC([O])C([O])C=C(22582)'],
    products = ['C#CC1OC[CH]C1[O](25774)'],
    transitionState = 'TS32',
    kinetics = Arrhenius(A=(2.07178e+07,'s^-1'), n=1.01592, Ea=(56.0365,'kJ/mol'), T0=(1,'K'), comment="""Estimated using template [R5_SS_D;doublebond_intra_pri_2H;radadd_intra] for rate rule [R5_SS_D;doublebond_intra_pri_2H;radadd_intra_O]
Euclidian distance = 1.0
family: Intra_R_Add_Endocyclic"""),
)

reaction(
    label = 'reaction33',
    reactants = ['C#CC([O])C([O])C=C(22582)'],
    products = ['C=CC1OC=[C]C1[O](25541)'],
    transitionState = 'TS33',
    kinetics = Arrhenius(A=(3.47e+11,'s^-1'), n=0.15, Ea=(58.576,'kJ/mol'), T0=(1,'K'), comment="""Estimated using template [R5_SS_T;triplebond_intra_H;radadd_intra] for rate rule [R5_SS_T;triplebond_intra_H;radadd_intra_O]
Euclidian distance = 1.0
family: Intra_R_Add_Endocyclic"""),
)

reaction(
    label = 'reaction34',
    reactants = ['C#CC([O])C([O])C=C(22582)'],
    products = ['C#CC(O)C(=O)C=C(25939)'],
    transitionState = 'TS34',
    kinetics = Arrhenius(A=(7.437e+08,'s^-1'), n=1.045, Ea=(63.4002,'kJ/mol'), T0=(1,'K'), comment="""Estimated using an average for rate rule [R3radExo;Y_rad;XH_Rrad]
Euclidian distance = 0
family: Intra_Disproportionation"""),
)

reaction(
    label = 'reaction35',
    reactants = ['C#CC([O])C([O])C=C(22582)'],
    products = ['C#CC(=O)C(O)C=C(25940)'],
    transitionState = 'TS35',
    kinetics = Arrhenius(A=(7.437e+08,'s^-1'), n=1.045, Ea=(63.4002,'kJ/mol'), T0=(1,'K'), comment="""Estimated using an average for rate rule [R3radExo;Y_rad;XH_Rrad]
Euclidian distance = 0
family: Intra_Disproportionation"""),
)

reaction(
    label = 'reaction36',
    reactants = ['C#CC([O])[C]([O])C[CH2](25941)'],
    products = ['C#CC([O])C([O])C=C(22582)'],
    transitionState = 'TS36',
    kinetics = Arrhenius(A=(4.48818e+10,'s^-1'), n=0.34095, Ea=(22.3009,'kJ/mol'), T0=(1,'K'), comment="""Estimated using average of templates [Rn;Y_rad_NDe;XH_Rrad] + [R2radExo;Y_rad;XH_Rrad] for rate rule [R2radExo;Y_rad_NDe;XH_Rrad]
Euclidian distance = 1.0
Multiplied by reaction path degeneracy 2.0
family: Intra_Disproportionation"""),
)

reaction(
    label = 'reaction37',
    reactants = ['[CH]=C=C([O])C([O])C[CH2](25942)'],
    products = ['C#CC([O])C([O])C=C(22582)'],
    transitionState = 'TS37',
    kinetics = Arrhenius(A=(5.2748e+09,'s^-1'), n=0.37, Ea=(88.9686,'kJ/mol'), T0=(1,'K'), comment="""Estimated using average of templates [R3;Y_rad_De;XH_Rrad] + [R3radExo;Y_rad;XH_Rrad] for rate rule [R3radExo;Y_rad_De;XH_Rrad]
Euclidian distance = 1.0
Multiplied by reaction path degeneracy 2.0
family: Intra_Disproportionation"""),
)

reaction(
    label = 'reaction38',
    reactants = ['[CH]=C=C([O])C([O])[CH]C(25943)'],
    products = ['C#CC([O])C([O])C=C(22582)'],
    transitionState = 'TS38',
    kinetics = Arrhenius(A=(1.34494e+09,'s^-1'), n=0.311, Ea=(39.225,'kJ/mol'), T0=(1,'K'), comment="""Estimated using template [R4;Y_rad_De;XH_Rrad] for rate rule [R4radEndo;Y_rad_De;XH_Rrad]
Euclidian distance = 1.0
Multiplied by reaction path degeneracy 3.0
family: Intra_Disproportionation"""),
)

reaction(
    label = 'reaction39',
    reactants = ['[C]#CC([O])C([O])C[CH2](25944)'],
    products = ['C#CC([O])C([O])C=C(22582)'],
    transitionState = 'TS39',
    kinetics = Arrhenius(A=(4.25221e+09,'s^-1'), n=0.137, Ea=(24.9733,'kJ/mol'), T0=(1,'K'), comment="""Estimated using template [R5;Y_rad;XH_Rrad] for rate rule [R5radExo;Y_rad;XH_Rrad]
Euclidian distance = 1.0
Multiplied by reaction path degeneracy 2.0
family: Intra_Disproportionation"""),
)

reaction(
    label = 'reaction40',
    reactants = ['[C]#CC([O])C([O])[CH]C(25945)'],
    products = ['C#CC([O])C([O])C=C(22582)'],
    transitionState = 'TS40',
    kinetics = Arrhenius(A=(6.37831e+09,'s^-1'), n=0.137, Ea=(24.9733,'kJ/mol'), T0=(1,'K'), comment="""Estimated using template [R6;Y_rad;XH_Rrad] for rate rule [R6radEndo;Y_rad;XH_Rrad]
Euclidian distance = 1.0
Multiplied by reaction path degeneracy 3.0
family: Intra_Disproportionation"""),
)

reaction(
    label = 'reaction41',
    reactants = ['C#CC([O])C([O])C=C(22582)'],
    products = ['C#CC1OOC1C=C(22598)'],
    transitionState = 'TS41',
    kinetics = Arrhenius(A=(1.62e+12,'s^-1'), n=-0.305, Ea=(8.28432,'kJ/mol'), T0=(1,'K'), comment="""Estimated using template [R4_SSS;Y_rad_out;Ypri_rad_out] for rate rule [R4_SSS;O_rad;Opri_rad]
Euclidian distance = 1.41421356237
family: Birad_recombination"""),
)

reaction(
    label = 'reaction46',
    reactants = ['C#CC([O])C([O])C=C(22582)'],
    products = ['[O]C=C=CCC=C[O](22588)'],
    transitionState = 'TS42',
    kinetics = Arrhenius(A=(2.214e+09,'s^-1'), n=0.749, Ea=(200.242,'kJ/mol'), T0=(1,'K'), comment="""From training reaction 2 used for hex_1_ene_5_yne
Exact match found for rate rule [hex_1_ene_5_yne]
Euclidian distance = 0
family: 6_membered_central_C-C_shift"""),
)

reaction(
    label = 'reaction43',
    reactants = ['[O]C1[C]=CC[CH]C1[O](25946)'],
    products = ['C#CC([O])C([O])C=C(22582)'],
    transitionState = 'TS43',
    kinetics = Arrhenius(A=(1e+13,'s^-1'), n=0, Ea=(0,'kJ/mol'), T0=(1,'K'), Tmin=(300,'K'), Tmax=(1500,'K'), comment="""Estimated using template [RJJ] for rate rule [R6JJ]
Euclidian distance = 1.0
family: 1,4_Cyclic_birad_scission"""),
)

reaction(
    label = 'reaction44',
    reactants = ['O(T)(63)', 'C#CC([O])[CH]C=C(23684)'],
    products = ['C#CC([O])C([O])C=C(22582)'],
    transitionState = 'TS44',
    kinetics = Arrhenius(A=(93609.6,'m^3/(mol*s)'), n=1.13083, Ea=(163.847,'kJ/mol'), T0=(1,'K'), Tmin=(303.03,'K'), Tmax=(2000,'K'), comment="""Estimated using template [Y_rad;O_birad] for rate rule [C_rad/H/OneDeC;O_birad]
Euclidian distance = 4.0
family: Birad_R_Recombination"""),
)

reaction(
    label = 'reaction45',
    reactants = ['O(T)(63)', '[CH]=C=CC([O])C=C(20057)'],
    products = ['C#CC([O])C([O])C=C(22582)'],
    transitionState = 'TS45',
    kinetics = Arrhenius(A=(93609.6,'m^3/(mol*s)'), n=1.13083, Ea=(163.847,'kJ/mol'), T0=(1,'K'), Tmin=(303.03,'K'), Tmax=(2000,'K'), comment="""Estimated using template [Y_rad;O_birad] for rate rule [C_rad/H/OneDeC;O_birad]
Euclidian distance = 4.0
family: Birad_R_Recombination"""),
)

network(
    label = '4883',
    isomers = [
        'C#CC([O])C([O])C=C(22582)',
    ],
    reactants = [
        ('C=CC=O(5269)', 'C#CC=O(21959)'),
    ],
    bathGas = {
        'N2': 0.25,
        'Ne': 0.25,
        'He': 0.25,
        'Ar': 0.25,
    },
)

pressureDependence(
    label = '4883',
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

