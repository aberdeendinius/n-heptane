species(
    label = 'C#CC([O])C([CH2])C=O(22594)',
    structure = SMILES('C#CC([O])C([CH2])C=O'),
    E0 = (281.373,'kJ/mol'),
    modes = [
        HarmonicOscillator(frequencies=([2782.5,750,1395,475,1775,1000,750,770,3400,2100,1380,1383.33,1386.67,1390,370,373.333,376.667,380,2800,3000,430,440,2175,525,3000,3100,440,815,1455,1000,388.962,388.962],'cm^-1')),
        HinderedRotor(inertia=(0.0974447,'amu*angstrom^2'), symmetry=1, barrier=(10.4616,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.0974439,'amu*angstrom^2'), symmetry=1, barrier=(10.4616,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.197874,'amu*angstrom^2'), symmetry=1, barrier=(21.2436,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.674186,'amu*angstrom^2'), symmetry=1, barrier=(72.3796,'kJ/mol'), semiclassical=False),
    ],
    spinMultiplicity = 3,
    opticalIsomers = 1,
    molecularWeight = (110.111,'amu'),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[0.337006,0.0819137,-9.87928e-05,6.33783e-08,-1.6199e-11,33972.1,30.6582], Tmin=(100,'K'), Tmax=(955.812,'K')), NASAPolynomial(coeffs=[13.7455,0.0257999,-1.07305e-05,1.95583e-09,-1.33398e-13,31408.9,-33.4241], Tmin=(955.812,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(281.373,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(315.95,'J/(mol*K)'), comment="""Thermo group additivity estimation: group(O2s-CsH) + group(Cs-(Cds-O2d)CsCsH) + group(Cs-(Cds-Cds)CsOsH) + group(Cs-CsHHH) + group(Cds-OdCsH) + group(Ct-CtCs) + group(Ct-CtH) + radical(CC(C)OJ) + radical(CJC(C)C=O)"""),
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
    label = 'C#CC([O])C1CC1[O](25458)',
    structure = SMILES('C#CC([O])C1CC1[O]'),
    E0 = (370.867,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
    molecularWeight = (110.111,'amu'),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[0.77377,0.0543621,-1.35851e-06,-4.85898e-08,2.58694e-11,44736.4,27.5852], Tmin=(100,'K'), Tmax=(935.864,'K')), NASAPolynomial(coeffs=[20.1807,0.0133854,-2.9523e-06,4.66399e-10,-3.63426e-14,39266,-74.5757], Tmin=(935.864,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(370.867,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(324.264,'J/(mol*K)'), comment="""Thermo group additivity estimation: group(O2s-CsH) + group(O2s-CsH) + group(Cs-CsCsCsH) + group(Cs-CsCsOsH) + group(Cs-CsCsHH) + group(Cs-(Cds-Cds)CsOsH) + group(Ct-CtCs) + group(Ct-CtH) + ring(Cyclopropane) + radical(CC(C)OJ) + radical(CC(C)OJ)"""),
)

species(
    label = '[CH]=C1OC1C([CH2])C=O(25459)',
    structure = SMILES('[CH]=C1OC1C([CH2])C=O'),
    E0 = (303.455,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
    molecularWeight = (110.111,'amu'),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[0.144777,0.0722331,-5.17441e-05,2.31562e-09,7.78978e-12,36647.4,28.1785], Tmin=(100,'K'), Tmax=(960.072,'K')), NASAPolynomial(coeffs=[21.6105,0.0120934,-3.55117e-06,6.32085e-10,-4.75902e-14,31175.6,-81.5376], Tmin=(960.072,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(303.455,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(320.107,'J/(mol*K)'), comment="""Thermo group additivity estimation: group(O2s-Cs(Cds-Cd)) + group(Cs-(Cds-O2d)CsCsH) + group(Cs-(Cds-Cds)CsOsH) + group(Cs-CsHHH) + group(Cds-CdsCsOs) + group(Cds-OdCsH) + group(Cds-CdsHH) + ring(methyleneoxirane) + radical(Cds_P) + radical(CJC(C)C=O)"""),
)

species(
    label = 'C#CC1OC([O])C1[CH2](25460)',
    structure = SMILES('C#CC1OC([O])C1[CH2]'),
    E0 = (334.608,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
    molecularWeight = (110.111,'amu'),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[0.0632304,0.0643852,-5.38673e-05,2.46679e-08,-4.28566e-12,40405,31.1408], Tmin=(100,'K'), Tmax=(1702.94,'K')), NASAPolynomial(coeffs=[12.6403,0.0210847,-3.60816e-06,2.48194e-10,-4.23941e-15,38116.4,-30.3743], Tmin=(1702.94,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(334.608,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(324.264,'J/(mol*K)'), comment="""Thermo group additivity estimation: group(O2s-CsCs) + group(O2s-CsH) + group(Cs-CsCsCsH) + group(Cs-(Cds-Cds)CsOsH) + group(Cs-CsOsOsH) + group(Cs-CsHHH) + group(Ct-CtCs) + group(Ct-CtH) + ring(Oxetane) + radical(Isobutyl) + radical(CCOJ)"""),
)

species(
    label = '[CH]=C1CC(C=O)C1[O](25461)',
    structure = SMILES('[CH]=C1CC(C=O)C1[O]'),
    E0 = (292.664,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
    molecularWeight = (110.111,'amu'),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[11.2057,-0.0119039,0.000120865,-1.17887e-07,2.6987e-11,34811.5,-16.5397], Tmin=(100,'K'), Tmax=(1745.61,'K')), NASAPolynomial(coeffs=[81.1828,0.018034,-6.83751e-05,1.68341e-08,-1.25082e-12,-18610.7,-476.165], Tmin=(1745.61,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(292.664,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(328.422,'J/(mol*K)'), comment="""Thermo group additivity estimation: group(O2s-CsH) + group(Cs-(Cds-O2d)CsCsH) + group(Cs-(Cds-Cds)CsOsH) + group(Cs-(Cds-Cds)CsHH) + group(Cds-CdsCsCs) + group(Cds-OdCsH) + group(Cds-CdsHH) + ring(methylenecyclobutane) + radical(CC(C)OJ) + radical(Cds_P)"""),
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
    label = 'C#CC([O])C(=C)C=O(25462)',
    structure = SMILES('C#CC([O])C(=C)C=O'),
    E0 = (176.422,'kJ/mol'),
    modes = [
        HarmonicOscillator(frequencies=([750,770,3400,2100,350,440,435,1725,2950,3100,1380,975,1025,1650,1380,1390,370,380,2900,435,2782.5,750,1395,475,1775,1000,2175,525,307.168,307.168],'cm^-1')),
        HinderedRotor(inertia=(0.338048,'amu*angstrom^2'), symmetry=1, barrier=(22.6364,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.338046,'amu*angstrom^2'), symmetry=1, barrier=(22.6364,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.338082,'amu*angstrom^2'), symmetry=1, barrier=(22.6363,'kJ/mol'), semiclassical=False),
    ],
    spinMultiplicity = 2,
    opticalIsomers = 1,
    molecularWeight = (109.103,'amu'),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[1.03811,0.067126,-6.58393e-05,3.37101e-08,-7.07509e-12,21323.6,26.5877], Tmin=(100,'K'), Tmax=(1129.72,'K')), NASAPolynomial(coeffs=[12.264,0.0273785,-1.30639e-05,2.5664e-09,-1.83156e-13,18787.2,-28.9397], Tmin=(1129.72,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(176.422,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(295.164,'J/(mol*K)'), comment="""Thermo group additivity estimation: group(O2s-CsH) + group(Cs-(Cds-Cds)(Cds-Cds)OsH) + group(Cd-CdCs(CO)) + group(Cds-CdsHH) + group(Cds-O2d(Cds-Cds)H) + group(Ct-CtCs) + group(Ct-CtH) + radical(CC(C)OJ)"""),
)

species(
    label = 'C#CC(=O)C([CH2])C=O(25463)',
    structure = SMILES('C#CC(=O)C([CH2])C=O'),
    E0 = (112.193,'kJ/mol'),
    modes = [
        HarmonicOscillator(frequencies=([750,770,3400,2100,3000,3100,440,815,1455,1000,375,552.5,462.5,1710,2782.5,750,1395,475,1775,1000,2175,525,1380,1390,370,380,2900,435,180],'cm^-1')),
        HinderedRotor(inertia=(0.00244726,'amu*angstrom^2'), symmetry=1, barrier=(7.94222,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(1.9265,'amu*angstrom^2'), symmetry=1, barrier=(44.2941,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(1.92294,'amu*angstrom^2'), symmetry=1, barrier=(44.2122,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.0369285,'amu*angstrom^2'), symmetry=1, barrier=(13.6106,'kJ/mol'), semiclassical=False),
    ],
    spinMultiplicity = 2,
    opticalIsomers = 1,
    molecularWeight = (109.103,'amu'),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[1.03107,0.0702204,-8.64097e-05,6.23052e-08,-1.87784e-11,13596.1,27.2665], Tmin=(100,'K'), Tmax=(799.16,'K')), NASAPolynomial(coeffs=[8.60396,0.0323134,-1.52541e-05,2.94221e-09,-2.06569e-13,12385.8,-7.56991], Tmin=(799.16,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(112.193,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(291.007,'J/(mol*K)'), comment="""Thermo group additivity estimation: group(Cs-CsCsCsH) + group(Cs-CsHHH) + group(Cds-O2d(Cds-Cds)Cs) + group(Cds-OdCsH) + group(Ct-CtCs) + group(Ct-CtH) + radical(CJC(C)C=O)"""),
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
    label = 'C#CC([O])C=C(13772)',
    structure = SMILES('C#CC([O])C=C'),
    E0 = (307.227,'kJ/mol'),
    modes = [
        HarmonicOscillator(frequencies=([2950,3100,1380,975,1025,1650,750,770,3400,2100,1380,1390,370,380,2900,435,2175,525,3010,987.5,1337.5,450,1655,217.186,217.21],'cm^-1')),
        HinderedRotor(inertia=(0.885817,'amu*angstrom^2'), symmetry=1, barrier=(29.6435,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.885441,'amu*angstrom^2'), symmetry=1, barrier=(29.6432,'kJ/mol'), semiclassical=False),
    ],
    spinMultiplicity = 2,
    opticalIsomers = 1,
    molecularWeight = (81.0926,'amu'),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[1.67511,0.043956,-2.37695e-05,-3.86807e-09,5.69169e-12,37041.1,21.6699], Tmin=(100,'K'), Tmax=(967.656,'K')), NASAPolynomial(coeffs=[12.4618,0.0157303,-5.38113e-06,9.38759e-10,-6.50844e-14,34187.4,-33.9737], Tmin=(967.656,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(307.227,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(249.434,'J/(mol*K)'), comment="""Thermo group additivity estimation: group(O2s-CsH) + group(Cs-(Cds-Cds)(Cds-Cds)OsH) + group(Cds-CdsCsH) + group(Cds-CdsHH) + group(Ct-CtCs) + group(Ct-CtH) + radical(CC(C)OJ)"""),
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
    label = '[CH2]C(C=O)C=O(12779)',
    structure = SMILES('[CH2]C(C=O)C=O'),
    E0 = (-104.853,'kJ/mol'),
    modes = [
        HarmonicOscillator(frequencies=([2695,2870,700,800,1380,1410,450,500,1750,1800,900,1100,1380,1390,370,380,2900,435,3000,3100,440,815,1455,1000],'cm^-1')),
        HinderedRotor(inertia=(1.36304,'amu*angstrom^2'), symmetry=1, barrier=(31.339,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.000540802,'amu*angstrom^2'), symmetry=1, barrier=(4.85933,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(1.36291,'amu*angstrom^2'), symmetry=1, barrier=(31.336,'kJ/mol'), semiclassical=False),
    ],
    spinMultiplicity = 2,
    opticalIsomers = 1,
    molecularWeight = (85.0813,'amu'),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[2.09608,0.0407753,-2.78866e-05,9.59164e-09,-1.36393e-12,-12541.8,21.6355], Tmin=(100,'K'), Tmax=(1580.53,'K')), NASAPolynomial(coeffs=[9.87979,0.0210762,-9.19104e-06,1.7058e-09,-1.16576e-13,-15002.3,-19.4794], Tmin=(1580.53,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(-104.853,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(245.277,'J/(mol*K)'), comment="""Thermo group additivity estimation: group(Cs-CsCsCsH) + group(Cs-CsHHH) + group(Cds-OdCsH) + group(Cds-OdCsH) + radical(CJC(C)C=O)"""),
)

species(
    label = 'C#CC([O])C(C)=C[O](25464)',
    structure = SMILES('C#CC([O])C(C)=C[O]'),
    E0 = (200.842,'kJ/mol'),
    modes = [
        HarmonicOscillator(frequencies=([750,770,3400,2100,2750,2800,2850,1350,1500,750,1050,1375,1000,350,440,435,1725,2175,525,1380,1390,370,380,2900,435,3010,987.5,1337.5,450,1655,180,180,180],'cm^-1')),
        HinderedRotor(inertia=(1.1896,'amu*angstrom^2'), symmetry=1, barrier=(27.3513,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(1.19095,'amu*angstrom^2'), symmetry=1, barrier=(27.3823,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(1.19093,'amu*angstrom^2'), symmetry=1, barrier=(27.3819,'kJ/mol'), semiclassical=False),
    ],
    spinMultiplicity = 3,
    opticalIsomers = 1,
    molecularWeight = (110.111,'amu'),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[-0.411889,0.0798789,-7.96809e-05,3.8838e-08,-7.17673e-12,24329.2,31.1966], Tmin=(100,'K'), Tmax=(1476.63,'K')), NASAPolynomial(coeffs=[21.473,0.0123654,-2.73875e-06,3.25788e-10,-1.74064e-14,18763.3,-79.8773], Tmin=(1476.63,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(200.842,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(320.107,'J/(mol*K)'), comment="""Thermo group additivity estimation: group(O2s-CsH) + group(O2s-(Cds-Cd)H) + group(Cs-(Cds-Cds)(Cds-Cds)OsH) + group(Cs-(Cds-Cds)HHH) + group(Cds-CdsCsCs) + group(Cds-CdsOsH) + group(Ct-CtCs) + group(Ct-CtH) + radical(C=COJ) + radical(CC(C)OJ)"""),
)

species(
    label = '[CH]=C=C(O)C([CH2])C=O(25465)',
    structure = SMILES('[CH]=C=C(O)C([CH2])C=O'),
    E0 = (161.011,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
    molecularWeight = (110.111,'amu'),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[-0.1913,0.0932431,-0.00012671,8.80238e-08,-2.38078e-11,19515.1,30.4896], Tmin=(100,'K'), Tmax=(913.504,'K')), NASAPolynomial(coeffs=[16.3191,0.0209404,-7.97447e-06,1.36197e-09,-8.84172e-14,16498.9,-47.6681], Tmin=(913.504,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(161.011,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(315.95,'J/(mol*K)'), comment="""Thermo group additivity estimation: group(O2s-(Cds-Cd)H) + group(Cs-(Cds-O2d)CsCsH) + group(Cs-CsHHH) + group(Cds-CdsCsOs) + group(Cds-OdCsH) + group(Cds-CdsHH) + group(Cdd-CdsCds) + radical(CJC(C)C=O) + radical(C=C=CJ)"""),
)

species(
    label = '[CH]=C=C([O])C(C)C=O(25466)',
    structure = SMILES('[CH]=C=C([O])C(C)C=O'),
    E0 = (88.3049,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
    molecularWeight = (110.111,'amu'),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[0.490686,0.0788846,-9.71378e-05,6.5646e-08,-1.77585e-11,10745.5,29.1937], Tmin=(100,'K'), Tmax=(903.885,'K')), NASAPolynomial(coeffs=[12.2619,0.026795,-1.0698e-05,1.89406e-09,-1.2636e-13,8617.43,-26.4064], Tmin=(903.885,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(88.3049,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(320.107,'J/(mol*K)'), comment="""Thermo group additivity estimation: group(O2s-(Cds-Cd)H) + group(Cs-(Cds-O2d)CsCsH) + group(Cs-CsHHH) + group(Cds-CdsCsOs) + group(Cds-OdCsH) + group(Cds-CdsHH) + group(Cdd-CdsCds) + radical(C=C(C)OJ) + radical(C=C=CJ)"""),
)

species(
    label = 'C#CC([O])C(C)[C]=O(25467)',
    structure = SMILES('C#CC([O])C(C)[C]=O'),
    E0 = (229.551,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
    molecularWeight = (110.111,'amu'),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[0.465738,0.0731921,-7.45331e-05,3.95856e-08,-8.33055e-12,27739.7,30.995], Tmin=(100,'K'), Tmax=(1157.9,'K')), NASAPolynomial(coeffs=[15.2652,0.0220669,-8.30317e-06,1.45346e-09,-9.75423e-14,24312.4,-42.5738], Tmin=(1157.9,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(229.551,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(315.95,'J/(mol*K)'), comment="""Thermo group additivity estimation: group(O2s-CsH) + group(Cs-(Cds-O2d)CsCsH) + group(Cs-(Cds-Cds)CsOsH) + group(Cs-CsHHH) + group(Cds-OdCsH) + group(Ct-CtCs) + group(Ct-CtH) + radical(CC(C)OJ) + radical(CC(C)CJ=O)"""),
)

species(
    label = 'C#CC(O)C([CH2])=C[O](25468)',
    structure = SMILES('C#CC(O)C([CH2])=C[O]'),
    E0 = (121.981,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
    molecularWeight = (110.111,'amu'),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[-0.808633,0.0840895,-8.63366e-05,4.23687e-08,-7.75282e-12,14862.6,33.0959], Tmin=(100,'K'), Tmax=(1540.91,'K')), NASAPolynomial(coeffs=[23.1767,0.00899053,-7.36056e-07,-7.19547e-11,9.91666e-15,8994.69,-88.0451], Tmin=(1540.91,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(121.981,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(315.95,'J/(mol*K)'), comment="""Thermo group additivity estimation: group(O2s-CsH) + group(O2s-(Cds-Cd)H) + group(Cs-(Cds-Cds)(Cds-Cds)OsH) + group(Cs-(Cds-Cds)HHH) + group(Cds-CdsCsCs) + group(Cds-CdsOsH) + group(Ct-CtCs) + group(Ct-CtH) + radical(Allyl_P) + radical(C=COJ)"""),
)

species(
    label = 'C#CC(O)C([CH2])[C]=O(25469)',
    structure = SMILES('C#CC(O)C([CH2])[C]=O'),
    E0 = (209.701,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
    molecularWeight = (110.111,'amu'),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[0.169102,0.0853422,-0.00010941,7.32258e-08,-1.92839e-11,25358.1,32.6074], Tmin=(100,'K'), Tmax=(933.567,'K')), NASAPolynomial(coeffs=[14.7363,0.0229285,-9.12958e-06,1.61612e-09,-1.07959e-13,22638.2,-36.6697], Tmin=(933.567,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(209.701,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(311.793,'J/(mol*K)'), comment="""Thermo group additivity estimation: group(O2s-CsH) + group(Cs-(Cds-O2d)CsCsH) + group(Cs-(Cds-Cds)CsOsH) + group(Cs-CsHHH) + group(Cds-OdCsH) + group(Ct-CtCs) + group(Ct-CtH) + radical(CC(C)CJ=O) + radical(CJC(C)C=O)"""),
)

species(
    label = '[C]#CC(O)C([CH2])C=O(25470)',
    structure = SMILES('[C]#CC(O)C([CH2])C=O'),
    E0 = (388.156,'kJ/mol'),
    modes = [
        HarmonicOscillator(frequencies=([3000,3100,440,815,1455,1000,2782.5,750,1395,475,1775,1000,3615,1277.5,1000,2175,525,1380,1383.33,1386.67,1390,370,373.333,376.667,380,2800,3000,430,440,280.054,280.059],'cm^-1')),
        HinderedRotor(inertia=(0.130264,'amu*angstrom^2'), symmetry=1, barrier=(7.24987,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.130275,'amu*angstrom^2'), symmetry=1, barrier=(7.24982,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.00214941,'amu*angstrom^2'), symmetry=1, barrier=(0.119627,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.248968,'amu*angstrom^2'), symmetry=1, barrier=(13.8567,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(1.16363,'amu*angstrom^2'), symmetry=1, barrier=(64.7699,'kJ/mol'), semiclassical=False),
    ],
    spinMultiplicity = 3,
    opticalIsomers = 1,
    molecularWeight = (110.111,'amu'),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[0.0596937,0.0911519,-0.000133608,1.04076e-07,-3.1597e-11,46822,32.582], Tmin=(100,'K'), Tmax=(902.004,'K')), NASAPolynomial(coeffs=[12.5848,0.0263921,-1.05889e-05,1.82528e-09,-1.17532e-13,44937.3,-24.4744], Tmin=(902.004,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(388.156,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(311.793,'J/(mol*K)'), comment="""Thermo group additivity estimation: group(O2s-CsH) + group(Cs-(Cds-O2d)CsCsH) + group(Cs-(Cds-Cds)CsOsH) + group(Cs-CsHHH) + group(Cds-OdCsH) + group(Ct-CtCs) + group(Ct-CtH) + radical(Acetyl) + radical(CJC(C)C=O)"""),
)

species(
    label = '[C]#CC([O])C(C)C=O(25471)',
    structure = SMILES('[C]#CC([O])C(C)C=O'),
    E0 = (408.006,'kJ/mol'),
    modes = [
        HarmonicOscillator(frequencies=([2750,2800,2850,1350,1500,750,1050,1375,1000,1380,1383.33,1386.67,1390,370,373.333,376.667,380,2800,3000,430,440,2175,525,2782.5,750,1395,475,1775,1000,332.719,340.751,353.559],'cm^-1')),
        HinderedRotor(inertia=(0.163751,'amu*angstrom^2'), symmetry=1, barrier=(13.1708,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.088778,'amu*angstrom^2'), symmetry=1, barrier=(7.6329,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.0014695,'amu*angstrom^2'), symmetry=1, barrier=(0.119627,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.832413,'amu*angstrom^2'), symmetry=1, barrier=(66.3863,'kJ/mol'), semiclassical=False),
    ],
    spinMultiplicity = 3,
    opticalIsomers = 1,
    molecularWeight = (110.111,'amu'),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[0.567956,0.0764015,-8.91886e-05,5.74996e-08,-1.49059e-11,49194.5,30.2179], Tmin=(100,'K'), Tmax=(940.575,'K')), NASAPolynomial(coeffs=[12.1361,0.0272063,-1.07346e-05,1.89332e-09,-1.26319e-13,47018.3,-24.8831], Tmin=(940.575,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(408.006,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(315.95,'J/(mol*K)'), comment="""Thermo group additivity estimation: group(O2s-CsH) + group(Cs-(Cds-O2d)CsCsH) + group(Cs-(Cds-Cds)CsOsH) + group(Cs-CsHHH) + group(Cds-OdCsH) + group(Ct-CtCs) + group(Ct-CtH) + radical(CC(C)OJ) + radical(Acetyl)"""),
)

species(
    label = '[CH]=[C]C([O])C=C(13781)',
    structure = SMILES('[CH]=[C]C([O])C=C'),
    E0 = (626.215,'kJ/mol'),
    modes = [
        HarmonicOscillator(frequencies=([1685,370,2950,3100,1380,975,1025,1650,1380,1390,370,380,2900,435,3120,650,792.5,1650,3010,987.5,1337.5,450,1655,373.097,374.84],'cm^-1')),
        HinderedRotor(inertia=(0.143307,'amu*angstrom^2'), symmetry=1, barrier=(14.1956,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.144824,'amu*angstrom^2'), symmetry=1, barrier=(14.2107,'kJ/mol'), semiclassical=False),
    ],
    spinMultiplicity = 4,
    opticalIsomers = 1,
    molecularWeight = (81.0926,'amu'),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[1.73288,0.0475138,-4.00564e-05,1.72941e-08,-3.02981e-12,75399.7,24.7241], Tmin=(100,'K'), Tmax=(1346.51,'K')), NASAPolynomial(coeffs=[11.3355,0.0189878,-8.27847e-06,1.56051e-09,-1.08621e-13,72813.8,-24.4597], Tmin=(1346.51,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(626.215,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(249.434,'J/(mol*K)'), comment="""Thermo group additivity estimation: group(O2s-CsH) + group(Cs-(Cds-Cds)(Cds-Cds)OsH) + group(Cds-CdsCsH) + group(Cds-CdsCsH) + group(Cds-CdsHH) + group(Cds-CdsHH) + radical(Cds_P) + radical(Cds_S) + radical(CC(C)OJ)"""),
)

species(
    label = 'C#CC([O])C([CH2])=C[O](25472)',
    structure = SMILES('C#CC([O])C([CH2])=C[O]'),
    E0 = (352.342,'kJ/mol'),
    modes = [
        HarmonicOscillator(frequencies=([750,770,3400,2100,350,440,435,1725,2175,525,3000,3100,440,815,1455,1000,3010,987.5,1337.5,450,1655,1380,1390,370,380,2900,435,180,180,180],'cm^-1')),
        HinderedRotor(inertia=(1.5717,'amu*angstrom^2'), symmetry=1, barrier=(36.1365,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(1.57248,'amu*angstrom^2'), symmetry=1, barrier=(36.1545,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(1.57244,'amu*angstrom^2'), symmetry=1, barrier=(36.1536,'kJ/mol'), semiclassical=False),
    ],
    spinMultiplicity = 4,
    opticalIsomers = 1,
    molecularWeight = (109.103,'amu'),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[-0.715804,0.0808591,-8.28394e-05,4.01315e-08,-7.21588e-12,42566.3,32.5752], Tmin=(100,'K'), Tmax=(1580.52,'K')), NASAPolynomial(coeffs=[22.9616,0.00737135,-2.21885e-07,-1.46863e-10,1.41048e-14,36776,-87.1331], Tmin=(1580.52,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(352.342,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(295.164,'J/(mol*K)'), comment="""Thermo group additivity estimation: group(O2s-CsH) + group(O2s-(Cds-Cd)H) + group(Cs-(Cds-Cds)(Cds-Cds)OsH) + group(Cs-(Cds-Cds)HHH) + group(Cds-CdsCsCs) + group(Cds-CdsOsH) + group(Ct-CtCs) + group(Ct-CtH) + radical(Allyl_P) + radical(C=COJ) + radical(CC(C)OJ)"""),
)

species(
    label = '[CH2]C([CH][O])C=O(13547)',
    structure = SMILES('[CH2]C([CH][O])C=O'),
    E0 = (223.822,'kJ/mol'),
    modes = [
        HarmonicOscillator(frequencies=([3025,407.5,1350,352.5,2782.5,750,1395,475,1775,1000,1380,1390,370,380,2900,435,3000,3100,440,815,1455,1000,180,1656.02],'cm^-1')),
        HinderedRotor(inertia=(0.236434,'amu*angstrom^2'), symmetry=1, barrier=(5.43608,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.237351,'amu*angstrom^2'), symmetry=1, barrier=(5.45718,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.238528,'amu*angstrom^2'), symmetry=1, barrier=(5.48423,'kJ/mol'), semiclassical=False),
    ],
    spinMultiplicity = 4,
    opticalIsomers = 1,
    molecularWeight = (85.0813,'amu'),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[1.2747,0.0721842,-0.000133189,1.30951e-07,-4.80698e-11,27006,24.7433], Tmin=(100,'K'), Tmax=(867.505,'K')), NASAPolynomial(coeffs=[3.3884,0.0336322,-1.67212e-05,3.1697e-09,-2.1486e-13,27723.2,21.0936], Tmin=(867.505,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(223.822,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(245.277,'J/(mol*K)'), comment="""Thermo group additivity estimation: group(O2s-CsH) + group(Cs-(Cds-O2d)CsCsH) + group(Cs-CsOsHH) + group(Cs-CsHHH) + group(Cds-OdCsH) + radical(CJC(C)C=O) + radical(CCsJOH) + radical(CCOJ)"""),
)

species(
    label = '[CH]=C=C([O])C([CH2])C=O(25473)',
    structure = SMILES('[CH]=C=C([O])C([CH2])C=O'),
    E0 = (298.816,'kJ/mol'),
    modes = [
        HarmonicOscillator(frequencies=([3120,650,792.5,1650,540,610,2055,350,440,435,1725,3000,3100,440,815,1455,1000,2782.5,750,1395,475,1775,1000,1380,1390,370,380,2900,435,180],'cm^-1')),
        HinderedRotor(inertia=(0.492968,'amu*angstrom^2'), symmetry=1, barrier=(11.3343,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.493133,'amu*angstrom^2'), symmetry=1, barrier=(11.3381,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.492969,'amu*angstrom^2'), symmetry=1, barrier=(11.3343,'kJ/mol'), semiclassical=False),
    ],
    spinMultiplicity = 4,
    opticalIsomers = 1,
    molecularWeight = (109.103,'amu'),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[0.145861,0.0896415,-0.000135691,1.07296e-07,-3.29169e-11,36073.4,30.7788], Tmin=(100,'K'), Tmax=(895.857,'K')), NASAPolynomial(coeffs=[12.6528,0.0241732,-9.95686e-06,1.7361e-09,-1.12375e-13,34218.7,-26.0288], Tmin=(895.857,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(298.816,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(295.164,'J/(mol*K)'), comment="""Thermo group additivity estimation: group(O2s-(Cds-Cd)H) + group(Cs-(Cds-O2d)CsCsH) + group(Cs-CsHHH) + group(Cds-CdsCsOs) + group(Cds-OdCsH) + group(Cds-CdsHH) + group(Cdd-CdsCds) + radical(CJC(C)C=O) + radical(C=C(C)OJ) + radical(C=C=CJ)"""),
)

species(
    label = 'C#CC([O])C([CH2])[C]=O(25474)',
    structure = SMILES('C#CC([O])C([CH2])[C]=O'),
    E0 = (440.061,'kJ/mol'),
    modes = [
        HarmonicOscillator(frequencies=([1855,455,950,750,770,3400,2100,1380,1383.33,1386.67,1390,370,373.333,376.667,380,2800,3000,430,440,2175,525,3000,3100,440,815,1455,1000,283.191,283.689],'cm^-1')),
        HinderedRotor(inertia=(0.00209412,'amu*angstrom^2'), symmetry=1, barrier=(0.119627,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.309148,'amu*angstrom^2'), symmetry=1, barrier=(17.6622,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.310587,'amu*angstrom^2'), symmetry=1, barrier=(17.6661,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(1.42818,'amu*angstrom^2'), symmetry=1, barrier=(81.5387,'kJ/mol'), semiclassical=False),
    ],
    spinMultiplicity = 4,
    opticalIsomers = 1,
    molecularWeight = (109.103,'amu'),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[0.341113,0.0812389,-0.000103116,6.76791e-08,-1.74557e-11,53058.3,31.7982], Tmin=(100,'K'), Tmax=(953.005,'K')), NASAPolynomial(coeffs=[14.6534,0.0211655,-8.56088e-06,1.53266e-09,-1.03306e-13,50330.4,-36.5612], Tmin=(953.005,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(440.061,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(291.007,'J/(mol*K)'), comment="""Thermo group additivity estimation: group(O2s-CsH) + group(Cs-(Cds-O2d)CsCsH) + group(Cs-(Cds-Cds)CsOsH) + group(Cs-CsHHH) + group(Cds-OdCsH) + group(Ct-CtCs) + group(Ct-CtH) + radical(CJC(C)C=O) + radical(CC(C)OJ) + radical(CC(C)CJ=O)"""),
)

species(
    label = '[C]#CC([O])C([CH2])C=O(25475)',
    structure = SMILES('[C]#CC([O])C([CH2])C=O'),
    E0 = (618.516,'kJ/mol'),
    modes = [
        HarmonicOscillator(frequencies=([2782.5,750,1395,475,1775,1000,1380,1383.33,1386.67,1390,370,373.333,376.667,380,2800,3000,430,440,2175,525,3000,3100,440,815,1455,1000,296.329,297.205,2795.26],'cm^-1')),
        HinderedRotor(inertia=(0.193967,'amu*angstrom^2'), symmetry=1, barrier=(12.2567,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.197921,'amu*angstrom^2'), symmetry=1, barrier=(12.2619,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(1.15663,'amu*angstrom^2'), symmetry=1, barrier=(72.3567,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(1.15545,'amu*angstrom^2'), symmetry=1, barrier=(72.3508,'kJ/mol'), semiclassical=False),
    ],
    spinMultiplicity = 4,
    opticalIsomers = 1,
    molecularWeight = (109.103,'amu'),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[0.225262,0.0871243,-0.000127581,9.88867e-08,-2.993e-11,74522.4,31.796], Tmin=(100,'K'), Tmax=(893.924,'K')), NASAPolynomial(coeffs=[12.4983,0.0246367,-1.00251e-05,1.7431e-09,-1.12993e-13,72630.6,-24.3465], Tmin=(893.924,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(618.516,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(291.007,'J/(mol*K)'), comment="""Thermo group additivity estimation: group(O2s-CsH) + group(Cs-(Cds-O2d)CsCsH) + group(Cs-(Cds-Cds)CsOsH) + group(Cs-CsHHH) + group(Cds-OdCsH) + group(Ct-CtCs) + group(Ct-CtH) + radical(Acetyl) + radical(CC(C)OJ) + radical(CJC(C)C=O)"""),
)

species(
    label = 'C#CC([O])C1[CH]OC1(25476)',
    structure = SMILES('C#CC([O])C1[CH]OC1'),
    E0 = (355.041,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
    molecularWeight = (110.111,'amu'),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[-0.857888,0.0774912,-7.5684e-05,3.73115e-08,-6.72557e-12,42901.8,31.9854], Tmin=(100,'K'), Tmax=(1675.94,'K')), NASAPolynomial(coeffs=[16.8877,0.0138477,3.3278e-07,-5.06586e-10,4.63979e-14,39943.6,-53.8702], Tmin=(1675.94,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(355.041,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(324.264,'J/(mol*K)'), comment="""Thermo group additivity estimation: group(O2s-CsCs) + group(O2s-CsH) + group(Cs-CsCsCsH) + group(Cs-(Cds-Cds)CsOsH) + group(Cs-CsOsHH) + group(Cs-CsOsHH) + group(Ct-CtCs) + group(Ct-CtH) + ring(Oxetane) + radical(CC(C)OJ) + radical(CCsJOCs)"""),
)

species(
    label = '[CH2]C(C=O)C1[C]=CO1(25388)',
    structure = SMILES('[CH2]C(C=O)C1[C]=CO1'),
    E0 = (274.735,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
    molecularWeight = (110.111,'amu'),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[0.327852,0.0640731,-2.23273e-05,-3.09224e-08,2.01143e-11,33190.4,28.1649], Tmin=(100,'K'), Tmax=(950.593,'K')), NASAPolynomial(coeffs=[22.7336,0.0106734,-2.57403e-06,4.65896e-10,-3.88262e-14,27083.6,-88.5102], Tmin=(950.593,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(274.735,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(320.107,'J/(mol*K)'), comment="""Thermo group additivity estimation: group(O2s-Cs(Cds-Cd)) + group(Cs-(Cds-O2d)CsCsH) + group(Cs-(Cds-Cds)CsOsH) + group(Cs-CsHHH) + group(Cds-CdsCsH) + group(Cds-OdCsH) + group(Cds-CdsOsH) + ring(Cyclobutene) + radical(Cds_S) + radical(CJC(C)C=O)"""),
)

species(
    label = 'C#CC1OO[CH]C1[CH2](25477)',
    structure = SMILES('C#CC1OO[CH]C1[CH2]'),
    E0 = (467.473,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
    molecularWeight = (110.111,'amu'),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[1.04411,0.0491009,1.79811e-05,-7.80688e-08,4.16265e-11,56345.9,25.3711], Tmin=(100,'K'), Tmax=(856.128,'K')), NASAPolynomial(coeffs=[19.6393,0.010011,2.73672e-06,-9.95327e-10,8.00809e-14,51410.5,-71.6803], Tmin=(856.128,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(467.473,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(324.264,'J/(mol*K)'), comment="""Thermo group additivity estimation: group(O2s-OsCs) + group(O2s-OsCs) + group(Cs-CsCsCsH) + group(Cs-(Cds-Cds)CsOsH) + group(Cs-CsOsHH) + group(Cs-CsHHH) + group(Ct-CtCs) + group(Ct-CtH) + ring(12dioxolane) + radical(CCsJOOC) + radical(Isobutyl)"""),
)

species(
    label = '[O]C1[C]=CCC1C=O(25478)',
    structure = SMILES('[O]C1[C]=CCC1C=O'),
    E0 = (224.718,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
    molecularWeight = (110.111,'amu'),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[1.39804,0.0438242,6.78485e-06,-3.90307e-08,1.73519e-11,27132.7,26.3071], Tmin=(100,'K'), Tmax=(1014.2,'K')), NASAPolynomial(coeffs=[14.2224,0.024136,-9.78381e-06,1.89201e-09,-1.38269e-13,22942.7,-43.5762], Tmin=(1014.2,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(224.718,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(328.422,'J/(mol*K)'), comment="""Thermo group additivity estimation: group(O2s-CsH) + group(Cs-(Cds-O2d)CsCsH) + group(Cs-(Cds-Cds)CsOsH) + group(Cs-(Cds-Cds)CsHH) + group(Cds-CdsCsH) + group(Cds-CdsCsH) + group(Cds-OdCsH) + ring(Cyclopentene) + radical(cyclopentene-vinyl) + radical(CC(C)OJ)"""),
)

species(
    label = 'C#CC(O)C(=C)C=O(25479)',
    structure = SMILES('C#CC(O)C(=C)C=O'),
    E0 = (-53.9392,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
    molecularWeight = (110.111,'amu'),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[0.882671,0.0710038,-7.12286e-05,3.79659e-08,-8.31908e-12,-6377.14,27.3397], Tmin=(100,'K'), Tmax=(1084.84,'K')), NASAPolynomial(coeffs=[12.1386,0.0295007,-1.3842e-05,2.69976e-09,-1.91974e-13,-8819.3,-27.8803], Tmin=(1084.84,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(-53.9392,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(315.95,'J/(mol*K)'), comment="""Thermo group additivity estimation: group(O2s-CsH) + group(Cs-(Cds-Cds)(Cds-Cds)OsH) + group(Cd-CdCs(CO)) + group(Cds-CdsHH) + group(Cds-O2d(Cds-Cds)H) + group(Ct-CtCs) + group(Ct-CtH)"""),
)

species(
    label = 'C#CC(=O)C(C)C=O(25480)',
    structure = SMILES('C#CC(=O)C(C)C=O'),
    E0 = (-98.3182,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
    molecularWeight = (110.111,'amu'),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[1.42562,0.059084,-4.73641e-05,2.08925e-08,-3.96844e-12,-11734.4,25.488], Tmin=(100,'K'), Tmax=(1196.37,'K')), NASAPolynomial(coeffs=[9.14549,0.0332732,-1.50029e-05,2.85963e-09,-2.00208e-13,-13581.6,-13.1401], Tmin=(1196.37,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(-98.3182,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(315.95,'J/(mol*K)'), comment="""Thermo group additivity estimation: group(Cs-CsCsCsH) + group(Cs-CsHHH) + group(Cds-O2d(Cds-Cds)Cs) + group(Cds-OdCsH) + group(Ct-CtCs) + group(Ct-CtH)"""),
)

species(
    label = 'C#CC([O])[C]([CH2])C[O](25481)',
    structure = SMILES('C#CC([O])[C]([CH2])C[O]'),
    E0 = (573.914,'kJ/mol'),
    modes = [
        HarmonicOscillator(frequencies=([360,370,350,750,770,3400,2100,2750,2850,1437.5,1250,1305,750,350,3000,3100,440,815,1455,1000,2175,525,1380,1390,370,380,2900,435,180,180,1033.05,1034.02],'cm^-1')),
        HinderedRotor(inertia=(0.00383524,'amu*angstrom^2'), symmetry=1, barrier=(2.90733,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.126473,'amu*angstrom^2'), symmetry=1, barrier=(2.90787,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(2.66015,'amu*angstrom^2'), symmetry=1, barrier=(61.1621,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.0807449,'amu*angstrom^2'), symmetry=1, barrier=(61.1212,'kJ/mol'), semiclassical=False),
    ],
    spinMultiplicity = 5,
    opticalIsomers = 1,
    molecularWeight = (110.111,'amu'),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[0.771328,0.0723589,-7.97282e-05,5.02068e-08,-1.29495e-11,69141,32.1143], Tmin=(100,'K'), Tmax=(937.975,'K')), NASAPolynomial(coeffs=[10.681,0.0300981,-1.21439e-05,2.17029e-09,-1.46052e-13,67282,-15.0596], Tmin=(937.975,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(573.914,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(315.95,'J/(mol*K)'), comment="""Thermo group additivity estimation: group(O2s-CsH) + group(O2s-CsH) + group(Cs-CsCsCsH) + group(Cs-(Cds-Cds)CsOsH) + group(Cs-CsOsHH) + group(Cs-CsHHH) + group(Ct-CtCs) + group(Ct-CtH) + radical(CCJ(C)CO) + radical(Isobutyl) + radical(CCOJ) + radical(CC(C)OJ)"""),
)

species(
    label = '[CH]=C=C([O])C([CH2])C[O](25482)',
    structure = SMILES('[CH]=C=C([O])C([CH2])C[O]'),
    E0 = (440.988,'kJ/mol'),
    modes = [
        HarmonicOscillator(frequencies=([3120,650,792.5,1650,540,610,2055,350,440,435,1725,2750,2850,1437.5,1250,1305,750,350,3000,3100,440,815,1455,1000,1380,1390,370,380,2900,435,180,180,915.491],'cm^-1')),
        HinderedRotor(inertia=(0.144146,'amu*angstrom^2'), symmetry=1, barrier=(3.3142,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.14098,'amu*angstrom^2'), symmetry=1, barrier=(3.24141,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.143101,'amu*angstrom^2'), symmetry=1, barrier=(3.29018,'kJ/mol'), semiclassical=False),
    ],
    spinMultiplicity = 5,
    opticalIsomers = 1,
    molecularWeight = (110.111,'amu'),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[0.372071,0.0827382,-0.000112066,8.37991e-08,-2.46604e-11,53166.4,33.4897], Tmin=(100,'K'), Tmax=(924.095,'K')), NASAPolynomial(coeffs=[11.4285,0.0274702,-1.0327e-05,1.72548e-09,-1.09302e-13,51439.3,-17.2667], Tmin=(924.095,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(440.988,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(320.107,'J/(mol*K)'), comment="""Thermo group additivity estimation: group(O2s-CsH) + group(O2s-(Cds-Cd)H) + group(Cs-(Cds-Cdd-Cd)CsCsH) + group(Cs-CsOsHH) + group(Cs-CsHHH) + group(Cds-CdsCsOs) + group(Cds-CdsHH) + group(Cdd-CdsCds) + radical(Isobutyl) + radical(CCOJ) + radical(C=C=CJ) + radical(C=C(C)OJ)"""),
)

species(
    label = '[CH]=C=C([O])C([CH2])[CH]O(25483)',
    structure = SMILES('[CH]=C=C([O])C([CH2])[CH]O'),
    E0 = (395.58,'kJ/mol'),
    modes = [
        HarmonicOscillator(frequencies=([3120,650,792.5,1650,3615,1277.5,1000,3025,407.5,1350,352.5,540,610,2055,350,440,435,1725,3000,3100,440,815,1455,1000,1380,1390,370,380,2900,435,313.979,316.634],'cm^-1')),
        HinderedRotor(inertia=(0.0939799,'amu*angstrom^2'), symmetry=1, barrier=(6.58793,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.0922529,'amu*angstrom^2'), symmetry=1, barrier=(6.6033,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.0922419,'amu*angstrom^2'), symmetry=1, barrier=(6.58576,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(1.31198,'amu*angstrom^2'), symmetry=1, barrier=(92.8304,'kJ/mol'), semiclassical=False),
    ],
    spinMultiplicity = 5,
    opticalIsomers = 1,
    molecularWeight = (110.111,'amu'),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[0.027256,0.0903171,-0.000114304,6.3699e-08,-8.74015e-12,47718,33.4469], Tmin=(100,'K'), Tmax=(720.888,'K')), NASAPolynomial(coeffs=[16.5849,0.0181503,-5.1475e-06,6.73754e-10,-3.42571e-14,44818.7,-44.5672], Tmin=(720.888,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(395.58,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(315.95,'J/(mol*K)'), comment="""Thermo group additivity estimation: group(O2s-CsH) + group(O2s-(Cds-Cd)H) + group(Cs-(Cds-Cdd-Cd)CsCsH) + group(Cs-CsOsHH) + group(Cs-CsHHH) + group(Cds-CdsCsOs) + group(Cds-CdsHH) + group(Cdd-CdsCds) + radical(C=C=CJ) + radical(Isobutyl) + radical(CCsJOH) + radical(C=C(C)OJ)"""),
)

species(
    label = '[C]#CC([O])C([CH2])C[O](25484)',
    structure = SMILES('[C]#CC([O])C([CH2])C[O]'),
    E0 = (758.484,'kJ/mol'),
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
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[0.406903,0.0839839,-0.000118982,9.41862e-08,-2.91978e-11,91349.2,33.9644], Tmin=(100,'K'), Tmax=(914.271,'K')), NASAPolynomial(coeffs=[9.88725,0.0306628,-1.20682e-05,2.05792e-09,-1.31604e-13,90110.7,-8.21592], Tmin=(914.271,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(758.484,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(315.95,'J/(mol*K)'), comment="""Thermo group additivity estimation: group(O2s-CsH) + group(O2s-CsH) + group(Cs-CsCsCsH) + group(Cs-(Cds-Cds)CsOsH) + group(Cs-CsOsHH) + group(Cs-CsHHH) + group(Ct-CtCs) + group(Ct-CtH) + radical(CC(C)OJ) + radical(Isobutyl) + radical(Acetyl) + radical(CCOJ)"""),
)

species(
    label = '[C]#CC([O])C([CH2])[CH]O(25485)',
    structure = SMILES('[C]#CC([O])C([CH2])[CH]O'),
    E0 = (713.077,'kJ/mol'),
    modes = [
        HarmonicOscillator(frequencies=([3000,3100,440,815,1455,1000,3025,407.5,1350,352.5,3615,1277.5,1000,2175,525,1380,1383.33,1386.67,1390,370,373.333,376.667,380,2800,3000,430,440,200,800,1200,1600],'cm^-1')),
        HinderedRotor(inertia=(0.156089,'amu*angstrom^2'), symmetry=1, barrier=(3.5888,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.156089,'amu*angstrom^2'), symmetry=1, barrier=(3.5888,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.156089,'amu*angstrom^2'), symmetry=1, barrier=(3.5888,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.156089,'amu*angstrom^2'), symmetry=1, barrier=(3.5888,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.156089,'amu*angstrom^2'), symmetry=1, barrier=(3.5888,'kJ/mol'), semiclassical=False),
    ],
    spinMultiplicity = 5,
    opticalIsomers = 1,
    molecularWeight = (110.111,'amu'),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[-0.242141,0.0957261,-0.000139038,1.03275e-07,-2.92832e-11,85913.8,34.9794], Tmin=(100,'K'), Tmax=(982.294,'K')), NASAPolynomial(coeffs=[15.0701,0.0212584,-6.82317e-06,9.87748e-10,-5.48424e-14,83490.1,-35.6445], Tmin=(982.294,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(713.077,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(311.793,'J/(mol*K)'), comment="""Thermo group additivity estimation: group(O2s-CsH) + group(O2s-CsH) + group(Cs-CsCsCsH) + group(Cs-(Cds-Cds)CsOsH) + group(Cs-CsOsHH) + group(Cs-CsHHH) + group(Ct-CtCs) + group(Ct-CtH) + radical(CC(C)OJ) + radical(Isobutyl) + radical(CCsJOH) + radical(Acetyl)"""),
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
    label = 'C#CC([O])C[CH2](22299)',
    structure = SMILES('C#CC([O])C[CH2]'),
    E0 = (388.448,'kJ/mol'),
    modes = [
        HarmonicOscillator(frequencies=([2750,2850,1437.5,1250,1305,750,350,750,770,3400,2100,1380,1390,370,380,2900,435,2175,525,3000,3100,440,815,1455,1000,537.043,942.673],'cm^-1')),
        HinderedRotor(inertia=(0.0258154,'amu*angstrom^2'), symmetry=1, barrier=(16.2745,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.707787,'amu*angstrom^2'), symmetry=1, barrier=(16.2734,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(3.69099,'amu*angstrom^2'), symmetry=1, barrier=(84.8631,'kJ/mol'), semiclassical=False),
    ],
    spinMultiplicity = 3,
    opticalIsomers = 1,
    molecularWeight = (82.1005,'amu'),
    collisionModel = TransportData(shapeIndex=2, epsilon=(3666,'J/mol'), sigma=(6.23254,'angstroms'), dipoleMoment=(0,'C*m'), polarizability=(0,'angstroms^3'), rotrelaxcollnum=0, comment="""Epsilon & sigma estimated with Tc=572.62 K, Pc=34.36 bar (from Joback method)"""),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[1.41119,0.0501189,-3.26148e-05,2.67964e-09,3.89556e-12,46818.9,24.0922], Tmin=(100,'K'), Tmax=(960.068,'K')), NASAPolynomial(coeffs=[13.001,0.0174105,-5.8526e-06,9.9829e-10,-6.79231e-14,43875.5,-35.0891], Tmin=(960.068,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(388.448,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(270.22,'J/(mol*K)'), comment="""Thermo group additivity estimation: group(O2s-CsH) + group(Cs-CsCsHH) + group(Cs-(Cds-Cds)CsOsH) + group(Cs-CsHHH) + group(Ct-CtCs) + group(Ct-CtH) + radical(RCCJ) + radical(CC(C)OJ)"""),
)

species(
    label = 'C#CC([O])CC=C[O](22586)',
    structure = SMILES('C#CC([O])CC=C[O]'),
    E0 = (219.874,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
    molecularWeight = (110.111,'amu'),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[0.349841,0.0658846,-3.00643e-05,-2.26308e-08,1.7727e-11,26589.5,30.1039], Tmin=(100,'K'), Tmax=(927.933,'K')), NASAPolynomial(coeffs=[21.5129,0.0116462,-2.17941e-06,2.92247e-10,-2.20858e-14,21069.4,-78.9933], Tmin=(927.933,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(219.874,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(320.107,'J/(mol*K)'), comment="""Thermo group additivity estimation: group(O2s-CsH) + group(O2s-(Cds-Cd)H) + group(Cs-(Cds-Cds)CsHH) + group(Cs-(Cds-Cds)CsOsH) + group(Cds-CdsCsH) + group(Cds-CdsOsH) + group(Ct-CtCs) + group(Ct-CtH) + radical(C=COJ) + radical(CC(C)OJ)"""),
)

species(
    label = 'C#CC([O])[CH]CC=O(25486)',
    structure = SMILES('C#CC([O])[CH]CC=O'),
    E0 = (276.525,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
    molecularWeight = (110.111,'amu'),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[0.723987,0.0698146,-6.86251e-05,3.5932e-08,-7.58446e-12,33378.1,31.547], Tmin=(100,'K'), Tmax=(1142.25,'K')), NASAPolynomial(coeffs=[13.4186,0.0253598,-1.0247e-05,1.85989e-09,-1.27199e-13,30478,-31.3854], Tmin=(1142.25,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(276.525,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(315.95,'J/(mol*K)'), comment="""Thermo group additivity estimation: group(O2s-CsH) + group(Cs-CsCsHH) + group(Cs-(Cds-Cds)CsOsH) + group(Cs-(Cds-O2d)CsHH) + group(Cds-OdCsH) + group(Ct-CtCs) + group(Ct-CtH) + radical(CCJCO) + radical(CC(C)OJ)"""),
)

species(
    label = 'C#CC1OCC1C=O(22601)',
    structure = SMILES('C#CC1OCC1C=O'),
    E0 = (24.533,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
    molecularWeight = (110.111,'amu'),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[1.27545,0.0487565,2.33776e-06,-4.88896e-08,2.7331e-11,3059.5,23.8479], Tmin=(100,'K'), Tmax=(862.399,'K')), NASAPolynomial(coeffs=[15.2155,0.017791,-2.40426e-06,7.71678e-11,3.80613e-15,-597.773,-48.6053], Tmin=(862.399,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(24.533,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(324.264,'J/(mol*K)'), comment="""Thermo group additivity estimation: group(O2s-CsCs) + group(Cs-(Cds-O2d)CsCsH) + group(Cs-(Cds-Cds)CsOsH) + group(Cs-CsOsHH) + group(Cds-OdCsH) + group(Ct-CtCs) + group(Ct-CtH) + ring(Oxetane)"""),
)

species(
    label = 'C#CC([O])OC=C[CH2](22590)',
    structure = SMILES('C#CC([O])OC=C[CH2]'),
    E0 = (224.258,'kJ/mol'),
    modes = [
        HarmonicOscillator(frequencies=([3000,3100,440,815,1455,1000,750,770,3400,2100,1380,1390,370,380,2900,435,2175,525,2995,3025,975,1000,1300,1375,400,500,1630,1680,269.297,269.306,269.311,269.312],'cm^-1')),
        HinderedRotor(inertia=(0.594159,'amu*angstrom^2'), symmetry=1, barrier=(30.5776,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.594174,'amu*angstrom^2'), symmetry=1, barrier=(30.5776,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.594101,'amu*angstrom^2'), symmetry=1, barrier=(30.5775,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.594155,'amu*angstrom^2'), symmetry=1, barrier=(30.5776,'kJ/mol'), semiclassical=False),
    ],
    spinMultiplicity = 3,
    opticalIsomers = 1,
    molecularWeight = (110.111,'amu'),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[-0.287086,0.0817232,-8.4109e-05,4.27363e-08,-8.34697e-12,27137,29.8495], Tmin=(100,'K'), Tmax=(1321.54,'K')), NASAPolynomial(coeffs=[20.8363,0.0146001,-4.30391e-06,6.527e-10,-4.06451e-14,21832.3,-76.895], Tmin=(1321.54,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(224.258,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(315.95,'J/(mol*K)'), comment="""Thermo group additivity estimation: group(O2s-Cs(Cds-Cd)) + group(O2s-CsH) + group(Cs-CsOsOsH) + group(Cs-(Cds-Cds)HHH) + group(Cds-CdsCsH) + group(Cds-CdsOsH) + group(Ct-CtCs) + group(Ct-CtH) + radical(Allyl_P) + radical(CCOJ)"""),
)

species(
    label = 'C#CC([O])C([CH2])=CO(25487)',
    structure = SMILES('C#CC([O])C([CH2])=CO'),
    E0 = (210.879,'kJ/mol'),
    modes = [
        HarmonicOscillator(frequencies=([750,770,3400,2100,3615,1277.5,1000,350,440,435,1725,3010,987.5,1337.5,450,1655,3000,3100,440,815,1455,1000,2175,525,1380,1390,370,380,2900,435,180,180],'cm^-1')),
        HinderedRotor(inertia=(1.30974,'amu*angstrom^2'), symmetry=1, barrier=(30.1136,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(1.30961,'amu*angstrom^2'), symmetry=1, barrier=(30.1105,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(1.3107,'amu*angstrom^2'), symmetry=1, barrier=(30.1356,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(1.31072,'amu*angstrom^2'), symmetry=1, barrier=(30.1361,'kJ/mol'), semiclassical=False),
    ],
    spinMultiplicity = 3,
    opticalIsomers = 1,
    molecularWeight = (110.111,'amu'),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[-1.3828,0.0908827,-9.71455e-05,4.84356e-08,-8.85188e-12,25580.3,33.5564], Tmin=(100,'K'), Tmax=(1585.63,'K')), NASAPolynomial(coeffs=[25.3103,0.00491196,1.80946e-06,-5.80353e-10,4.47397e-14,19457.7,-100.141], Tmin=(1585.63,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(210.879,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(315.95,'J/(mol*K)'), comment="""Thermo group additivity estimation: group(O2s-CsH) + group(O2s-(Cds-Cd)H) + group(Cs-(Cds-Cds)(Cds-Cds)OsH) + group(Cs-(Cds-Cds)HHH) + group(Cds-CdsCsCs) + group(Cds-CdsOsH) + group(Ct-CtCs) + group(Ct-CtH) + radical(CC(C)OJ) + radical(Allyl_P)"""),
)

species(
    label = '[CH2]C1[CH]OC=[C]C1[O](25488)',
    structure = SMILES('[CH2]C1[CH]OC=[C]C1[O]'),
    E0 = (507.197,'kJ/mol'),
    modes = [
        HarmonicOscillator(frequencies=([2750,2883.33,3016.67,3150,900,966.667,1033.33,1100,3000,3100,440,815,1455,1000,300,800,800,800,800,800,800,800,800,800,800,1600,1600,1600,1600,1600,1600,1600,1600,1600,1600],'cm^-1')),
        HinderedRotor(inertia=(0.156089,'amu*angstrom^2'), symmetry=1, barrier=(3.5888,'kJ/mol'), semiclassical=False),
    ],
    spinMultiplicity = 5,
    opticalIsomers = 1,
    molecularWeight = (110.111,'amu'),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[0.721502,0.0520527,1.63814e-05,-7.86078e-08,4.06117e-11,61138.7,25.3461], Tmin=(100,'K'), Tmax=(897.355,'K')), NASAPolynomial(coeffs=[23.3847,0.00559182,2.8409e-06,-7.91003e-10,5.53091e-14,54874.6,-93.7764], Tmin=(897.355,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(507.197,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(328.422,'J/(mol*K)'), comment="""Thermo group additivity estimation: group(O2s-Cs(Cds-Cd)) + group(O2s-CsH) + group(Cs-CsCsCsH) + group(Cs-(Cds-Cds)CsOsH) + group(Cs-CsOsHH) + group(Cs-CsHHH) + group(Cds-CdsCsH) + group(Cds-CdsOsH) + ring(3,4-Dihydro-2H-pyran) + radical(Isobutyl) + radical(CCsJOC(O)) + radical(CC(C)OJ) + radical(Cds_S)"""),
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
    label = 'C#CC([O])C=C[O](23769)',
    structure = SMILES('C#CC([O])C=C[O]'),
    E0 = (239.897,'kJ/mol'),
    modes = [
        HarmonicOscillator(frequencies=([750,770,3400,2100,1380,1390,370,380,2900,435,2175,525,2995,3025,975,1000,1300,1375,400,500,1630,1680,180,180,180],'cm^-1')),
        HinderedRotor(inertia=(1.48012,'amu*angstrom^2'), symmetry=1, barrier=(34.0308,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(1.48163,'amu*angstrom^2'), symmetry=1, barrier=(34.0656,'kJ/mol'), semiclassical=False),
    ],
    spinMultiplicity = 3,
    opticalIsomers = 1,
    molecularWeight = (96.0841,'amu'),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[1.06484,0.0525536,-2.4138e-05,-1.98503e-08,1.54148e-11,28969.8,25.1307], Tmin=(100,'K'), Tmax=(921.344,'K')), NASAPolynomial(coeffs=[18.874,0.00662152,-4.56314e-07,-1.20305e-11,-7.40412e-16,24356,-66.5587], Tmin=(921.344,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(239.897,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(249.434,'J/(mol*K)'), comment="""Thermo group additivity estimation: group(O2s-CsH) + group(O2s-(Cds-Cd)H) + group(Cs-(Cds-Cds)(Cds-Cds)OsH) + group(Cds-CdsCsH) + group(Cds-CdsOsH) + group(Ct-CtCs) + group(Ct-CtH) + radical(CC(C)OJ) + radical(C=COJ)"""),
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
    label = '[CH]=C=CC([CH2])C=O(19966)',
    structure = SMILES('[CH]=C=CC([CH2])C=O'),
    E0 = (375.569,'kJ/mol'),
    modes = [
        HarmonicOscillator(frequencies=([3120,650,792.5,1650,540,610,2055,3000,3100,440,815,1455,1000,2782.5,750,1395,475,1775,1000,3010,987.5,1337.5,450,1655,1380,1390,370,380,2900,435],'cm^-1')),
        HinderedRotor(inertia=(0.532983,'amu*angstrom^2'), symmetry=1, barrier=(12.2543,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.531765,'amu*angstrom^2'), symmetry=1, barrier=(12.2263,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(1.058,'amu*angstrom^2'), symmetry=1, barrier=(24.3256,'kJ/mol'), semiclassical=False),
    ],
    spinMultiplicity = 3,
    opticalIsomers = 1,
    molecularWeight = (94.1112,'amu'),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[0.835179,0.0742316,-0.000103291,8.24322e-08,-2.64427e-11,45280.1,26.7472], Tmin=(100,'K'), Tmax=(843.498,'K')), NASAPolynomial(coeffs=[8.98071,0.0294866,-1.2842e-05,2.34628e-09,-1.58011e-13,44123.6,-9.87394], Tmin=(843.498,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(375.569,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(295.164,'J/(mol*K)'), comment="""Thermo group additivity estimation: group(Cs-(Cds-O2d)CsCsH) + group(Cs-CsHHH) + group(Cds-CdsCsH) + group(Cds-OdCsH) + group(Cds-CdsHH) + group(Cdd-CdsCds) + radical(CJC(C)C=O) + radical(C=C=CJ)"""),
)

species(
    label = '[CH]C(C=O)C([O])C#C(25489)',
    structure = SMILES('[CH]C(C=O)C([O])C#C'),
    E0 = (519.077,'kJ/mol'),
    modes = [
        HarmonicOscillator(frequencies=([2782.5,750,1395,475,1775,1000,750,770,3400,2100,1380,1383.33,1386.67,1390,370,373.333,376.667,380,2800,3000,430,440,2175,525,200,800,1066.67,1333.33,1600],'cm^-1')),
        HinderedRotor(inertia=(0.156089,'amu*angstrom^2'), symmetry=1, barrier=(3.5888,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.156089,'amu*angstrom^2'), symmetry=1, barrier=(3.5888,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.156089,'amu*angstrom^2'), symmetry=1, barrier=(3.5888,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.156089,'amu*angstrom^2'), symmetry=1, barrier=(3.5888,'kJ/mol'), semiclassical=False),
    ],
    spinMultiplicity = 4,
    opticalIsomers = 1,
    molecularWeight = (109.103,'amu'),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[0.475824,0.0761377,-8.70131e-05,5.15432e-08,-1.20546e-11,62558.7,30.8086], Tmin=(100,'K'), Tmax=(1046.28,'K')), NASAPolynomial(coeffs=[14.7962,0.0213887,-8.52055e-06,1.52844e-09,-1.03777e-13,59562.1,-38.9267], Tmin=(1046.28,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(519.077,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(291.007,'J/(mol*K)'), comment="""Thermo group additivity estimation: group(O2s-CsH) + group(Cs-(Cds-O2d)CsCsH) + group(Cs-(Cds-Cds)CsOsH) + group(Cs-CsHHH) + group(Cds-OdCsH) + group(Ct-CtCs) + group(Ct-CtH) + radical(CC(C)OJ) + radical(CCJ2_triplet)"""),
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
    E0 = (281.373,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS2',
    E0 = (370.867,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS3',
    E0 = (365.818,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS4',
    E0 = (403.416,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS5',
    E0 = (353.032,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS6',
    E0 = (399.044,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS7',
    E0 = (356.632,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS8',
    E0 = (281.373,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS9',
    E0 = (360.25,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS10',
    E0 = (281.373,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS11',
    E0 = (505.307,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS12',
    E0 = (404.184,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS13',
    E0 = (435.281,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS14',
    E0 = (399.312,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS15',
    E0 = (399.312,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS16',
    E0 = (356.653,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS17',
    E0 = (408.766,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS18',
    E0 = (472.045,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS19',
    E0 = (447.905,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS20',
    E0 = (360.168,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS21',
    E0 = (659.577,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS22',
    E0 = (564.146,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS23',
    E0 = (775.758,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS24',
    E0 = (510.62,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS25',
    E0 = (651.866,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS26',
    E0 = (830.321,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS27',
    E0 = (467.798,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS28',
    E0 = (523.543,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS29',
    E0 = (467.473,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS30',
    E0 = (339.949,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS31',
    E0 = (344.773,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS32',
    E0 = (344.773,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS33',
    E0 = (596.215,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS34',
    E0 = (529.956,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS35',
    E0 = (434.805,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS36',
    E0 = (783.457,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS37',
    E0 = (738.05,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS38',
    E0 = (993.645,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS39',
    E0 = (441.308,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS40',
    E0 = (404.123,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS41',
    E0 = (289.657,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS42',
    E0 = (538.058,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS43',
    E0 = (354.746,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS44',
    E0 = (507.197,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS45',
    E0 = (655.583,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS46',
    E0 = (782.451,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS47',
    E0 = (730.882,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

reaction(
    label = 'reaction1',
    reactants = ['C#CC([O])C([CH2])C=O(22594)'],
    products = ['C=CC=O(5269)', 'C#CC=O(21959)'],
    transitionState = 'TS1',
    kinetics = Arrhenius(A=(5e+12,'s^-1'), n=0, Ea=(0,'kJ/mol'), T0=(1,'K'), Tmin=(300,'K'), Tmax=(1500,'K'), comment="""Exact match found for rate rule [RJJ]
Euclidian distance = 0
family: 1,4_Linear_birad_scission"""),
)

reaction(
    label = 'reaction2',
    reactants = ['C#CC([O])C([CH2])C=O(22594)'],
    products = ['C#CC([O])C1CC1[O](25458)'],
    transitionState = 'TS2',
    kinetics = Arrhenius(A=(3.93521e+09,'s^-1'), n=0.743095, Ea=(89.4945,'kJ/mol'), T0=(1,'K'), comment="""Estimated using template [R4_S;multiplebond_intra;radadd_intra_cs2H] for rate rule [R4_S_CO;carbonylbond_intra_H;radadd_intra_cs2H]
Euclidian distance = 2.2360679775
family: Intra_R_Add_Exocyclic
Ea raised from 85.9 to 89.5 kJ/mol to match endothermicity of reaction."""),
)

reaction(
    label = 'reaction3',
    reactants = ['C#CC([O])C([CH2])C=O(22594)'],
    products = ['[CH]=C1OC1C([CH2])C=O(25459)'],
    transitionState = 'TS3',
    kinetics = Arrhenius(A=(3.98513e+09,'s^-1'), n=0.768, Ea=(84.4457,'kJ/mol'), T0=(1,'K'), comment="""Estimated using template [R4_S_T;triplebond_intra_H;radadd_intra] for rate rule [R4_S_T;triplebond_intra_H;radadd_intra_O]
Euclidian distance = 1.0
family: Intra_R_Add_Exocyclic"""),
)

reaction(
    label = 'reaction4',
    reactants = ['C#CC([O])C([CH2])C=O(22594)'],
    products = ['C#CC1OC([O])C1[CH2](25460)'],
    transitionState = 'TS4',
    kinetics = Arrhenius(A=(2.724e+10,'s^-1'), n=0.478, Ea=(122.043,'kJ/mol'), T0=(1,'K'), comment="""Estimated using template [R5_SS;multiplebond_intra;radadd_intra_O] for rate rule [R5_SS_CO;carbonylbond_intra_H;radadd_intra_O]
Euclidian distance = 2.2360679775
family: Intra_R_Add_Exocyclic"""),
)

reaction(
    label = 'reaction5',
    reactants = ['C#CC([O])C([CH2])C=O(22594)'],
    products = ['[CH]=C1CC(C=O)C1[O](25461)'],
    transitionState = 'TS5',
    kinetics = Arrhenius(A=(177207,'s^-1'), n=1.88643, Ea=(71.6589,'kJ/mol'), T0=(1,'K'), comment="""Estimated using template [R5_SS;multiplebond_intra;radadd_intra_cs2H] for rate rule [R5_SS_T;triplebond_intra_H;radadd_intra_cs2H]
Euclidian distance = 2.2360679775
family: Intra_R_Add_Exocyclic"""),
)

reaction(
    label = 'reaction6',
    reactants = ['H(8)', 'C#CC([O])C(=C)C=O(25462)'],
    products = ['C#CC([O])C([CH2])C=O(22594)'],
    transitionState = 'TS6',
    kinetics = Arrhenius(A=(72.1434,'m^3/(mol*s)'), n=1.66666, Ea=(10.8177,'kJ/mol'), T0=(1,'K'), comment="""Estimated using template [Cds-OneDeCs_Cds;HJ] for rate rule [Cds-COCs_Cds;HJ]
Euclidian distance = 1.0
family: R_Addition_MultipleBond"""),
)

reaction(
    label = 'reaction7',
    reactants = ['H(8)', 'C#CC(=O)C([CH2])C=O(25463)'],
    products = ['C#CC([O])C([CH2])C=O(22594)'],
    transitionState = 'TS7',
    kinetics = Arrhenius(A=(1.59e+07,'cm^3/(mol*s)'), n=1.84, Ea=(32.6352,'kJ/mol'), T0=(1,'K'), Tmin=(300,'K'), Tmax=(1500,'K'), comment="""From training reaction 2830 used for CO-CtCs_O;HJ
Exact match found for rate rule [CO-CtCs_O;HJ]
Euclidian distance = 0
family: R_Addition_MultipleBond"""),
)

reaction(
    label = 'reaction8',
    reactants = ['C=CC=O(5269)', '[CH]=C=C[O](8556)'],
    products = ['C#CC([O])C([CH2])C=O(22594)'],
    transitionState = 'TS8',
    kinetics = Arrhenius(A=(1.76206,'m^3/(mol*s)'), n=1.97634, Ea=(92.8367,'kJ/mol'), T0=(1,'K'), comment="""Estimated using average of templates [Cds-OneDeH_Cds;CJ] + [Cds-COH_Cds;YJ] for rate rule [Cds-COH_Cds;CJ]
Euclidian distance = 1.0
family: R_Addition_MultipleBond
Ea raised from 89.8 to 92.8 kJ/mol to match endothermicity of reaction."""),
)

reaction(
    label = 'reaction9',
    reactants = ['[CH]=O(373)', 'C#CC([O])C=C(13772)'],
    products = ['C#CC([O])C([CH2])C=O(22594)'],
    transitionState = 'TS9',
    kinetics = Arrhenius(A=(0.00168615,'m^3/(mol*s)'), n=2.52599, Ea=(19.6608,'kJ/mol'), T0=(1,'K'), comment="""Estimated using template [Cds-CsH_Cds-HH;CJ] for rate rule [Cds-CsH_Cds-HH;CO_pri_rad]
Euclidian distance = 2.0
family: R_Addition_MultipleBond"""),
)

reaction(
    label = 'reaction10',
    reactants = ['[CH2]C=C[O](5266)', 'C#CC=O(21959)'],
    products = ['C#CC([O])C([CH2])C=O(22594)'],
    transitionState = 'TS10',
    kinetics = Arrhenius(A=(2.99e+06,'cm^3/(mol*s)'), n=2.12, Ea=(106.786,'kJ/mol'), T0=(1,'K'), Tmin=(300,'K'), Tmax=(1500,'K'), comment="""Estimated using template [CO-CtH_O;YJ] for rate rule [CO-CtH_O;CJ]
Euclidian distance = 1.0
family: R_Addition_MultipleBond
Ea raised from 103.8 to 106.8 kJ/mol to match endothermicity of reaction."""),
)

reaction(
    label = 'reaction11',
    reactants = ['[C]#C(5143)', '[CH2]C(C=O)C=O(12779)'],
    products = ['C#CC([O])C([CH2])C=O(22594)'],
    transitionState = 'TS11',
    kinetics = Arrhenius(A=(4679.9,'m^3/(mol*s)'), n=0.573452, Ea=(58.2237,'kJ/mol'), T0=(1,'K'), comment="""Estimated using template [CO_O;CJ] for rate rule [CO-CsH_O;CtJ_Ct]
Euclidian distance = 2.82842712475
Multiplied by reaction path degeneracy 2.0
family: R_Addition_MultipleBond"""),
)

reaction(
    label = 'reaction12',
    reactants = ['C#CC([O])C(C)=C[O](25464)'],
    products = ['C#CC([O])C([CH2])C=O(22594)'],
    transitionState = 'TS12',
    kinetics = Arrhenius(A=(2.307e+09,'s^-1'), n=1.31, Ea=(203.342,'kJ/mol'), T0=(1,'K'), Tmin=(300,'K'), Tmax=(1500,'K'), comment="""From training reaction 163 used for R2H_S;C_rad_out_OneDe/Cs;Cs_H_out_2H
Exact match found for rate rule [R2H_S;C_rad_out_OneDe/Cs;Cs_H_out_2H]
Euclidian distance = 0
Multiplied by reaction path degeneracy 3.0
family: intra_H_migration"""),
)

reaction(
    label = 'reaction13',
    reactants = ['C#CC([O])C([CH2])C=O(22594)'],
    products = ['[CH]=C=C(O)C([CH2])C=O(25465)'],
    transitionState = 'TS13',
    kinetics = Arrhenius(A=(1.70223e+09,'s^-1'), n=1.15155, Ea=(153.908,'kJ/mol'), T0=(1,'K'), comment="""Estimated using template [R2H_S;Y_rad_out;Cs_H_out_OneDe] for rate rule [R2H_S;O_rad_out;Cs_H_out_Ct]
Euclidian distance = 1.41421356237
family: intra_H_migration"""),
)

reaction(
    label = 'reaction14',
    reactants = ['C#CC([O])C([CH2])C=O(22594)'],
    products = ['[CH]=C=C([O])C(C)C=O(25466)'],
    transitionState = 'TS14',
    kinetics = Arrhenius(A=(83345.1,'s^-1'), n=2.17519, Ea=(117.939,'kJ/mol'), T0=(1,'K'), comment="""Estimated using an average for rate rule [R3H_SS_Cs;C_rad_out_2H;XH_out]
Euclidian distance = 0
family: intra_H_migration"""),
)

reaction(
    label = 'reaction15',
    reactants = ['C#CC([O])C([CH2])C=O(22594)'],
    products = ['C#CC([O])C(C)[C]=O(25467)'],
    transitionState = 'TS15',
    kinetics = Arrhenius(A=(83345.1,'s^-1'), n=2.17519, Ea=(117.939,'kJ/mol'), T0=(1,'K'), comment="""Estimated using template [R3H_SS_Cs;C_rad_out_2H;XH_out] for rate rule [R3H_SS_Cs;C_rad_out_2H;CO_H_out]
Euclidian distance = 1.0
family: intra_H_migration"""),
)

reaction(
    label = 'reaction16',
    reactants = ['C#CC([O])C([CH2])C=O(22594)'],
    products = ['C#CC(O)C([CH2])=C[O](25468)'],
    transitionState = 'TS16',
    kinetics = Arrhenius(A=(111914,'s^-1'), n=2.27675, Ea=(75.2806,'kJ/mol'), T0=(1,'K'), comment="""Estimated using an average for rate rule [R3H_SS_Cs;O_rad_out;XH_out]
Euclidian distance = 0
family: intra_H_migration"""),
)

reaction(
    label = 'reaction17',
    reactants = ['C#CC([O])C([CH2])C=O(22594)'],
    products = ['C#CC(O)C([CH2])[C]=O(25469)'],
    transitionState = 'TS17',
    kinetics = Arrhenius(A=(1.75172e+06,'s^-1'), n=1.80068, Ea=(127.394,'kJ/mol'), T0=(1,'K'), comment="""Estimated using template [R4H_SSS;O_rad_out;XH_out] for rate rule [R4H_SSS;O_rad_out;CO_H_out]
Euclidian distance = 1.0
family: intra_H_migration"""),
)

reaction(
    label = 'reaction18',
    reactants = ['[C]#CC(O)C([CH2])C=O(25470)'],
    products = ['C#CC([O])C([CH2])C=O(22594)'],
    transitionState = 'TS18',
    kinetics = Arrhenius(A=(0.00504539,'s^-1'), n=3.83, Ea=(83.8892,'kJ/mol'), T0=(1,'K'), comment="""Estimated using template [R4H_RSS;Y_rad_out;O_H_out] for rate rule [R4H_TSS;Ct_rad_out;O_H_out]
Euclidian distance = 1.41421356237
family: intra_H_migration"""),
)

reaction(
    label = 'reaction19',
    reactants = ['[C]#CC([O])C(C)C=O(25471)'],
    products = ['C#CC([O])C([CH2])C=O(22594)'],
    transitionState = 'TS19',
    kinetics = Arrhenius(A=(263079,'s^-1'), n=1.73643, Ea=(39.8993,'kJ/mol'), T0=(1,'K'), comment="""Estimated using template [R5H_RSSR;Y_rad_out;Cs_H_out_2H] for rate rule [R5H_TSSS;Ct_rad_out;Cs_H_out_2H]
Euclidian distance = 2.2360679775
Multiplied by reaction path degeneracy 3.0
family: intra_H_migration"""),
)

reaction(
    label = 'reaction20',
    reactants = ['[CH2]C=C[O](5266)', '[CH]=C=C[O](8556)'],
    products = ['C#CC([O])C([CH2])C=O(22594)'],
    transitionState = 'TS20',
    kinetics = Arrhenius(A=(1.9789e+07,'m^3/(mol*s)'), n=-0.126319, Ea=(0,'kJ/mol'), T0=(1,'K'), comment="""Estimated using an average for rate rule [Y_rad;Y_rad]
Euclidian distance = 0
family: R_Recombination
Ea raised from -15.6 to -15.6 kJ/mol.
Ea raised from -15.6 to 0 kJ/mol."""),
)

reaction(
    label = 'reaction21',
    reactants = ['[CH]=O(373)', '[CH]=[C]C([O])C=C(13781)'],
    products = ['C#CC([O])C([CH2])C=O(22594)'],
    transitionState = 'TS21',
    kinetics = Arrhenius(A=(3.7839e+06,'m^3/(mol*s)'), n=0, Ea=(0,'kJ/mol'), T0=(1,'K'), comment="""Estimated using an average for rate rule [Y_rad;CO_pri_rad]
Euclidian distance = 0
family: R_Recombination
Ea raised from -3.8 to 0 kJ/mol."""),
)

reaction(
    label = 'reaction22',
    reactants = ['H(8)', 'C#CC([O])C([CH2])=C[O](25472)'],
    products = ['C#CC([O])C([CH2])C=O(22594)'],
    transitionState = 'TS22',
    kinetics = Arrhenius(A=(4.34078e+06,'m^3/(mol*s)'), n=0.278577, Ea=(0,'kJ/mol'), T0=(1,'K'), comment="""Estimated using an average for rate rule [Y_rad;H_rad]
Euclidian distance = 0
family: R_Recombination
Ea raised from -1.4 to 0 kJ/mol."""),
)

reaction(
    label = 'reaction23',
    reactants = ['[C]#C(5143)', '[CH2]C([CH][O])C=O(13547)'],
    products = ['C#CC([O])C([CH2])C=O(22594)'],
    transitionState = 'TS23',
    kinetics = Arrhenius(A=(1.34536e+08,'m^3/(mol*s)'), n=0, Ea=(0,'kJ/mol'), T0=(1,'K'), comment="""Estimated using an average for rate rule [Ct_rad/Ct;Y_rad]
Euclidian distance = 0
family: R_Recombination"""),
)

reaction(
    label = 'reaction24',
    reactants = ['H(8)', '[CH]=C=C([O])C([CH2])C=O(25473)'],
    products = ['C#CC([O])C([CH2])C=O(22594)'],
    transitionState = 'TS24',
    kinetics = Arrhenius(A=(4.34078e+06,'m^3/(mol*s)'), n=0.278577, Ea=(0,'kJ/mol'), T0=(1,'K'), comment="""Estimated using an average for rate rule [Y_rad;H_rad]
Euclidian distance = 0
family: R_Recombination
Ea raised from -1.4 to 0 kJ/mol."""),
)

reaction(
    label = 'reaction25',
    reactants = ['H(8)', 'C#CC([O])C([CH2])[C]=O(25474)'],
    products = ['C#CC([O])C([CH2])C=O(22594)'],
    transitionState = 'TS25',
    kinetics = Arrhenius(A=(1.53107e+07,'m^3/(mol*s)'), n=-0.133333, Ea=(0,'kJ/mol'), T0=(1,'K'), comment="""Estimated using template [CO_rad/NonDe;Y_rad] for rate rule [CO_rad/NonDe;H_rad]
Euclidian distance = 1.0
family: R_Recombination"""),
)

reaction(
    label = 'reaction26',
    reactants = ['H(8)', '[C]#CC([O])C([CH2])C=O(25475)'],
    products = ['C#CC([O])C([CH2])C=O(22594)'],
    transitionState = 'TS26',
    kinetics = Arrhenius(A=(1.81e+14,'cm^3/(mol*s)','*|/',3), n=0, Ea=(0,'kJ/mol'), T0=(1,'K'), Tmin=(300,'K'), Tmax=(2500,'K'), comment="""From training reaction 61 used for H_rad;Ct_rad/Ct
Exact match found for rate rule [H_rad;Ct_rad/Ct]
Euclidian distance = 0
family: R_Recombination"""),
)

reaction(
    label = 'reaction27',
    reactants = ['C#CC([O])C([CH2])C=O(22594)'],
    products = ['C#CC([O])C1[CH]OC1(25476)'],
    transitionState = 'TS27',
    kinetics = Arrhenius(A=(7.01137e+09,'s^-1'), n=0.572544, Ea=(186.425,'kJ/mol'), T0=(1,'K'), comment="""Estimated using average of templates [R4_S;multiplebond_intra;radadd_intra_cs2H] + [R4_S_CO;carbonyl_intra;radadd_intra] for rate rule [R4_S_CO;carbonyl_intra_H;radadd_intra_cs2H]
Euclidian distance = 2.2360679775
family: Intra_R_Add_Endocyclic"""),
)

reaction(
    label = 'reaction28',
    reactants = ['C#CC([O])C([CH2])C=O(22594)'],
    products = ['[CH2]C(C=O)C1[C]=CO1(25388)'],
    transitionState = 'TS28',
    kinetics = Arrhenius(A=(1.503e+11,'s^-1'), n=0.221, Ea=(242.17,'kJ/mol'), T0=(1,'K'), comment="""Estimated using template [R4_S;multiplebond_intra;radadd_intra_O] for rate rule [R4_S_T;triplebond_intra_H;radadd_intra_O]
Euclidian distance = 2.2360679775
family: Intra_R_Add_Endocyclic"""),
)

reaction(
    label = 'reaction29',
    reactants = ['C#CC([O])C([CH2])C=O(22594)'],
    products = ['C#CC1OO[CH]C1[CH2](25477)'],
    transitionState = 'TS29',
    kinetics = Arrhenius(A=(2.81184e+09,'s^-1'), n=0.551229, Ea=(186.101,'kJ/mol'), T0=(1,'K'), comment="""Estimated using template [R5_SS;multiplebond_intra;radadd_intra] for rate rule [R5_SS_CO;carbonyl_intra_H;radadd_intra_O]
Euclidian distance = 2.44948974278
family: Intra_R_Add_Endocyclic
Ea raised from 182.1 to 186.1 kJ/mol to match endothermicity of reaction."""),
)

reaction(
    label = 'reaction30',
    reactants = ['C#CC([O])C([CH2])C=O(22594)'],
    products = ['[O]C1[C]=CCC1C=O(25478)'],
    transitionState = 'TS30',
    kinetics = Arrhenius(A=(3.47e+11,'s^-1'), n=0.15, Ea=(58.576,'kJ/mol'), T0=(1,'K'), comment="""Estimated using template [R5_SS_T;triplebond_intra_H;radadd_intra] for rate rule [R5_SS_T;triplebond_intra_H;radadd_intra_cs2H]
Euclidian distance = 2.0
family: Intra_R_Add_Endocyclic"""),
)

reaction(
    label = 'reaction31',
    reactants = ['C#CC([O])C([CH2])C=O(22594)'],
    products = ['C#CC(O)C(=C)C=O(25479)'],
    transitionState = 'TS31',
    kinetics = Arrhenius(A=(7.437e+08,'s^-1'), n=1.045, Ea=(63.4002,'kJ/mol'), T0=(1,'K'), comment="""Estimated using an average for rate rule [R3radExo;Y_rad;XH_Rrad]
Euclidian distance = 0
family: Intra_Disproportionation"""),
)

reaction(
    label = 'reaction32',
    reactants = ['C#CC([O])C([CH2])C=O(22594)'],
    products = ['C#CC(=O)C(C)C=O(25480)'],
    transitionState = 'TS32',
    kinetics = Arrhenius(A=(7.437e+08,'s^-1'), n=1.045, Ea=(63.4002,'kJ/mol'), T0=(1,'K'), comment="""Estimated using an average for rate rule [R3radExo;Y_rad;XH_Rrad]
Euclidian distance = 0
family: Intra_Disproportionation"""),
)

reaction(
    label = 'reaction33',
    reactants = ['C#CC([O])[C]([CH2])C[O](25481)'],
    products = ['C#CC([O])C([CH2])C=O(22594)'],
    transitionState = 'TS33',
    kinetics = Arrhenius(A=(4.48818e+10,'s^-1'), n=0.34095, Ea=(22.3009,'kJ/mol'), T0=(1,'K'), comment="""Estimated using average of templates [Rn;Y_rad_NDe;XH_Rrad] + [R2radExo;Y_rad;XH_Rrad] for rate rule [R2radExo;Y_rad_NDe;XH_Rrad]
Euclidian distance = 1.0
Multiplied by reaction path degeneracy 2.0
family: Intra_Disproportionation"""),
)

reaction(
    label = 'reaction34',
    reactants = ['[CH]=C=C([O])C([CH2])C[O](25482)'],
    products = ['C#CC([O])C([CH2])C=O(22594)'],
    transitionState = 'TS34',
    kinetics = Arrhenius(A=(5.2748e+09,'s^-1'), n=0.37, Ea=(88.9686,'kJ/mol'), T0=(1,'K'), comment="""Estimated using average of templates [R3;Y_rad_De;XH_Rrad] + [R3radExo;Y_rad;XH_Rrad] for rate rule [R3radExo;Y_rad_De;XH_Rrad]
Euclidian distance = 1.0
Multiplied by reaction path degeneracy 2.0
family: Intra_Disproportionation"""),
)

reaction(
    label = 'reaction35',
    reactants = ['[CH]=C=C([O])C([CH2])[CH]O(25483)'],
    products = ['C#CC([O])C([CH2])C=O(22594)'],
    transitionState = 'TS35',
    kinetics = Arrhenius(A=(4.48312e+08,'s^-1'), n=0.311, Ea=(39.225,'kJ/mol'), T0=(1,'K'), comment="""Estimated using template [R4;Y_rad_De;XH_Rrad] for rate rule [R4radEndo;Y_rad_De;XH_Rrad]
Euclidian distance = 1.0
family: Intra_Disproportionation"""),
)

reaction(
    label = 'reaction36',
    reactants = ['[C]#CC([O])C([CH2])C[O](25484)'],
    products = ['C#CC([O])C([CH2])C=O(22594)'],
    transitionState = 'TS36',
    kinetics = Arrhenius(A=(4.25221e+09,'s^-1'), n=0.137, Ea=(24.9733,'kJ/mol'), T0=(1,'K'), comment="""Estimated using template [R5;Y_rad;XH_Rrad] for rate rule [R5radExo;Y_rad;XH_Rrad]
Euclidian distance = 1.0
Multiplied by reaction path degeneracy 2.0
family: Intra_Disproportionation"""),
)

reaction(
    label = 'reaction37',
    reactants = ['[C]#CC([O])C([CH2])[CH]O(25485)'],
    products = ['C#CC([O])C([CH2])C=O(22594)'],
    transitionState = 'TS37',
    kinetics = Arrhenius(A=(2.1261e+09,'s^-1'), n=0.137, Ea=(24.9733,'kJ/mol'), T0=(1,'K'), comment="""Estimated using template [R6;Y_rad;XH_Rrad] for rate rule [R6radEndo;Y_rad;XH_Rrad]
Euclidian distance = 1.0
family: Intra_Disproportionation"""),
)

reaction(
    label = 'reaction38',
    reactants = ['[C-]#[O+](374)', 'C#CC([O])C[CH2](22299)'],
    products = ['C#CC([O])C([CH2])C=O(22594)'],
    transitionState = 'TS38',
    kinetics = Arrhenius(A=(0.118397,'m^3/(mol*s)'), n=2.3675, Ea=(305.306,'kJ/mol'), T0=(1,'K'), comment="""Estimated using an average for rate rule [CO;R_H]
Euclidian distance = 0
Multiplied by reaction path degeneracy 2.0
family: 1,2_Insertion_CO"""),
)

reaction(
    label = 'reaction39',
    reactants = ['C#CC([O])C([CH2])C=O(22594)'],
    products = ['C#CC([O])CC=C[O](22586)'],
    transitionState = 'TS39',
    kinetics = Arrhenius(A=(6.55606e+10,'s^-1'), n=0.64, Ea=(159.935,'kJ/mol'), T0=(1,'K'), comment="""Estimated using template [cCs(-HC)CJ;CsJ;C] for rate rule [cCs(-HC)CJ;CsJ-HH;C]
Euclidian distance = 1.0
family: 1,2_shiftC"""),
)

reaction(
    label = 'reaction40',
    reactants = ['C#CC([O])C([CH2])C=O(22594)'],
    products = ['C#CC([O])[CH]CC=O(25486)'],
    transitionState = 'TS40',
    kinetics = Arrhenius(A=(8.889e+11,'s^-1'), n=0.232, Ea=(122.75,'kJ/mol'), T0=(1,'K'), comment="""Estimated using template [cCs(-HC)CJ;CsJ;CO] for rate rule [cCs(-HC)CJ;CsJ-HH;CO]
Euclidian distance = 1.0
family: 1,2_shiftC"""),
)

reaction(
    label = 'reaction41',
    reactants = ['C#CC([O])C([CH2])C=O(22594)'],
    products = ['C#CC1OCC1C=O(22601)'],
    transitionState = 'TS41',
    kinetics = Arrhenius(A=(1.62e+12,'s^-1'), n=-0.305, Ea=(8.28432,'kJ/mol'), T0=(1,'K'), Tmin=(600,'K'), Tmax=(2000,'K'), comment="""Estimated using template [R4_SSS;C_rad_out_2H;Ypri_rad_out] for rate rule [R4_SSS;C_rad_out_2H;Opri_rad]
Euclidian distance = 1.0
family: Birad_recombination"""),
)

reaction(
    label = 'reaction42',
    reactants = ['C#CC([O])OC=C[CH2](22590)'],
    products = ['C#CC([O])C([CH2])C=O(22594)'],
    transitionState = 'TS42',
    kinetics = Arrhenius(A=(7040,'s^-1'), n=2.66, Ea=(313.8,'kJ/mol'), T0=(1,'K'), Tmin=(300,'K'), Tmax=(1500,'K'), comment="""Estimated using an average for rate rule [R_ROR;R1_doublebond;R2_doublebond_H;R_O_C]
Euclidian distance = 0
family: ketoenol"""),
)

reaction(
    label = 'reaction43',
    reactants = ['C#CC([O])C([CH2])=CO(25487)'],
    products = ['C#CC([O])C([CH2])C=O(22594)'],
    transitionState = 'TS43',
    kinetics = Arrhenius(A=(605.045,'s^-1'), n=2.96, Ea=(143.867,'kJ/mol'), T0=(1,'K'), comment="""Estimated using an average for rate rule [R_ROR;R1_doublebond;R2_doublebond_H;R_O_H]
Euclidian distance = 0
family: ketoenol"""),
)

reaction(
    label = 'reaction44',
    reactants = ['[CH2]C1[CH]OC=[C]C1[O](25488)'],
    products = ['C#CC([O])C([CH2])C=O(22594)'],
    transitionState = 'TS44',
    kinetics = Arrhenius(A=(1e+13,'s^-1'), n=0, Ea=(0,'kJ/mol'), T0=(1,'K'), Tmin=(300,'K'), Tmax=(1500,'K'), comment="""Estimated using template [RJJ] for rate rule [R6JJ]
Euclidian distance = 1.0
family: 1,4_Cyclic_birad_scission"""),
)

reaction(
    label = 'reaction45',
    reactants = ['CH2(T)(28)', 'C#CC([O])C=C[O](23769)'],
    products = ['C#CC([O])C([CH2])C=O(22594)'],
    transitionState = 'TS45',
    kinetics = Arrhenius(A=(1.14854e+06,'m^3/(mol*s)'), n=0.575199, Ea=(34.3157,'kJ/mol'), T0=(1,'K'), comment="""Estimated using template [Y_rad;Birad] for rate rule [C_rad/H/OneDeC;Birad]
Euclidian distance = 4.0
family: Birad_R_Recombination"""),
)

reaction(
    label = 'reaction46',
    reactants = ['O(T)(63)', '[CH]=C=CC([CH2])C=O(19966)'],
    products = ['C#CC([O])C([CH2])C=O(22594)'],
    transitionState = 'TS46',
    kinetics = Arrhenius(A=(93609.6,'m^3/(mol*s)'), n=1.13083, Ea=(163.847,'kJ/mol'), T0=(1,'K'), Tmin=(303.03,'K'), Tmax=(2000,'K'), comment="""Estimated using template [Y_rad;O_birad] for rate rule [C_rad/H/OneDeC;O_birad]
Euclidian distance = 4.0
family: Birad_R_Recombination"""),
)

reaction(
    label = 'reaction47',
    reactants = ['H(8)', '[CH]C(C=O)C([O])C#C(25489)'],
    products = ['C#CC([O])C([CH2])C=O(22594)'],
    transitionState = 'TS47',
    kinetics = Arrhenius(A=(1e+07,'m^3/(mol*s)'), n=0, Ea=(0,'kJ/mol'), T0=(1,'K'), comment="""Estimated using an average for rate rule [H_rad;Birad]
Euclidian distance = 0
family: Birad_R_Recombination"""),
)

network(
    label = '4895',
    isomers = [
        'C#CC([O])C([CH2])C=O(22594)',
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
    label = '4895',
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

