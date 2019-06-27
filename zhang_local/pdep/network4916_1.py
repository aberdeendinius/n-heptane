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
    label = '[CH]=[C]C1OC1[C]=C(25025)',
    structure = SMILES('[CH]=[C]C1OC1[C]=C'),
    E0 = (790.257,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
    molecularWeight = (93.1033,'amu'),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[0.330787,0.061028,-5.77598e-05,2.78653e-08,-5.00211e-12,95195.1,29.749], Tmin=(100,'K'), Tmax=(1638.66,'K')), NASAPolynomial(coeffs=[14.8222,0.0124228,-1.15563e-06,-9.08818e-11,1.47579e-14,92222.2,-41.9005], Tmin=(1638.66,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(790.257,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(274.378,'J/(mol*K)'), comment="""Thermo group additivity estimation: group(O2s-CsCs) + group(Cs-(Cds-Cds)CsOsH) + group(Cs-(Cds-Cds)CsOsH) + group(Cds-CdsCsH) + group(Cds-CdsCsH) + group(Cds-CdsHH) + group(Cds-CdsHH) + ring(Ethylene_oxide) + radical(Cds_S) + radical(Cds_S) + radical(Cds_P)"""),
)

species(
    label = '[CH]=[C]C1OC=C1[CH2](25026)',
    structure = SMILES('[CH]=[C]C1OC=C1[CH2]'),
    E0 = (660.13,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
    molecularWeight = (93.1033,'amu'),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[0.979164,0.0443145,2.44544e-05,-8.13215e-08,3.91019e-11,79524.5,23.0897], Tmin=(100,'K'), Tmax=(932.957,'K')), NASAPolynomial(coeffs=[24.1021,0.00169879,2.09552e-06,-4.06887e-10,1.85643e-14,72750.1,-100.044], Tmin=(932.957,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(660.13,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(274.378,'J/(mol*K)'), comment="""Thermo group additivity estimation: group(O2s-Cs(Cds-Cd)) + group(Cs-(Cds-Cds)(Cds-Cds)OsH) + group(Cs-(Cds-Cds)HHH) + group(Cds-CdsCsCs) + group(Cds-CdsCsH) + group(Cds-CdsOsH) + group(Cds-CdsHH) + ring(Cyclobutene) + radical(Cds_P) + radical(Cds_S) + radical(Allyl_P)"""),
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
    label = 'C=C=C[O](12570)',
    structure = SMILES('C=C=C[O]'),
    E0 = (115.398,'kJ/mol'),
    modes = [
        HarmonicOscillator(frequencies=([2950,3100,1380,975,1025,1650,3010,987.5,1337.5,450,1655,540,610,2055,180],'cm^-1')),
    ],
    spinMultiplicity = 2,
    opticalIsomers = 1,
    molecularWeight = (55.0553,'amu'),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[2.73497,0.0165704,2.41315e-05,-5.11524e-08,2.31238e-11,13935.2,12.0298], Tmin=(100,'K'), Tmax=(921.695,'K')), NASAPolynomial(coeffs=[12.9299,0.00148105,1.24089e-06,-2.76514e-10,1.55495e-14,10817.5,-43.0414], Tmin=(921.695,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(115.398,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(157.975,'J/(mol*K)'), comment="""Thermo group additivity estimation: group(O2s-(Cds-Cd)H) + group(Cds-CdsOsH) + group(Cds-CdsHH) + group(Cdd-CdsCds) + radical(C=COJ)"""),
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
    label = '[CH]=C=CO[CH]C=[CH](24816)',
    structure = SMILES('[CH]=C=CO[CH]C=[CH]'),
    E0 = (577.143,'kJ/mol'),
    modes = [
        HarmonicOscillator(frequencies=([3115,3125,620,680,785,800,1600,1700,3025,407.5,1350,352.5,2995,3025,975,1000,1300,1375,400,500,1630,1680,540,610,2055,180,180],'cm^-1')),
        HinderedRotor(inertia=(1.48331,'amu*angstrom^2'), symmetry=1, barrier=(34.1042,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(1.48289,'amu*angstrom^2'), symmetry=1, barrier=(34.0945,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(1.48434,'amu*angstrom^2'), symmetry=1, barrier=(34.1279,'kJ/mol'), semiclassical=False),
    ],
    spinMultiplicity = 4,
    opticalIsomers = 1,
    molecularWeight = (93.1033,'amu'),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[0.12046,0.0692171,-7.01415e-05,3.449e-08,-6.37828e-12,69567.7,28.5432], Tmin=(100,'K'), Tmax=(1503.81,'K')), NASAPolynomial(coeffs=[19.2299,0.00941307,-1.53706e-06,1.07847e-10,-2.6868e-15,64835.1,-68.0712], Tmin=(1503.81,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(577.143,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(270.22,'J/(mol*K)'), comment="""Thermo group additivity estimation: group(O2s-Cs(Cds-Cd)) + group(Cs-(Cds-Cds)OsHH) + group(Cds-CdsCsH) + group(Cds-CdsOsH) + group(Cds-CdsHH) + group(Cds-CdsHH) + group(Cdd-CdsCds) + radical(C=C=CJ) + radical(C=CCJ(O)C) + radical(Cds_P)"""),
)

species(
    label = '[CH]=C=[C]OC[C]=C(24326)',
    structure = SMILES('[CH]=C=[C]OC[C]=C'),
    E0 = (696.693,'kJ/mol'),
    modes = [
        HarmonicOscillator(frequencies=([1670,1700,300,440,2750,2850,1437.5,1250,1305,750,350,540,610,2055,3120,650,792.5,1650,2950,3100,1380,975,1025,1650,180,180,180],'cm^-1')),
        HinderedRotor(inertia=(0.860325,'amu*angstrom^2'), symmetry=1, barrier=(19.7806,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.860382,'amu*angstrom^2'), symmetry=1, barrier=(19.7819,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.860393,'amu*angstrom^2'), symmetry=1, barrier=(19.7821,'kJ/mol'), semiclassical=False),
    ],
    spinMultiplicity = 4,
    opticalIsomers = 1,
    molecularWeight = (93.1033,'amu'),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[0.754126,0.0675478,-7.68102e-05,4.49661e-08,-1.02651e-11,83913,28.5234], Tmin=(100,'K'), Tmax=(1079.73,'K')), NASAPolynomial(coeffs=[14.6325,0.0161354,-5.38891e-06,8.6948e-10,-5.53591e-14,80916,-39.4968], Tmin=(1079.73,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(696.693,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(270.22,'J/(mol*K)'), comment="""Thermo group additivity estimation: group(O2s-Cs(Cds-Cd)) + group(Cs-(Cds-Cds)OsHH) + group(Cds-CdsCsH) + group(Cds-CdsOsH) + group(Cds-CdsHH) + group(Cds-CdsHH) + group(Cdd-CdsCds) + radical(C=CJO) + radical(Cds_S) + radical(C=C=CJ)"""),
)

species(
    label = '[CH]=[C]COC=C=[CH](22625)',
    structure = SMILES('[CH]=[C]COC=C=[CH]'),
    E0 = (704.045,'kJ/mol'),
    modes = [
        HarmonicOscillator(frequencies=([3115,3125,620,680,785,800,1600,1700,2750,2850,1437.5,1250,1305,750,350,540,610,2055,3010,987.5,1337.5,450,1655,1685,370,180,180],'cm^-1')),
        HinderedRotor(inertia=(1.13356,'amu*angstrom^2'), symmetry=1, barrier=(26.0627,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(1.13434,'amu*angstrom^2'), symmetry=1, barrier=(26.0807,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(1.13417,'amu*angstrom^2'), symmetry=1, barrier=(26.0767,'kJ/mol'), semiclassical=False),
    ],
    spinMultiplicity = 4,
    opticalIsomers = 1,
    molecularWeight = (93.1033,'amu'),
    collisionModel = TransportData(shapeIndex=2, epsilon=(3451.63,'J/mol'), sigma=(5.81624,'angstroms'), dipoleMoment=(0,'C*m'), polarizability=(0,'angstroms^3'), rotrelaxcollnum=0, comment="""Epsilon & sigma estimated with Tc=539.14 K, Pc=39.81 bar (from Joback method)"""),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[0.0799838,0.0739297,-8.28071e-05,4.47525e-08,-9.08264e-12,84828.8,28.2162], Tmin=(100,'K'), Tmax=(1353.79,'K')), NASAPolynomial(coeffs=[19.5873,0.00814717,-8.9536e-07,-2.86214e-11,7.6414e-15,80293.4,-69.0473], Tmin=(1353.79,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(704.045,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(270.22,'J/(mol*K)'), comment="""Thermo group additivity estimation: group(O2s-Cs(Cds-Cd)) + group(Cs-(Cds-Cds)OsHH) + group(Cds-CdsCsH) + group(Cds-CdsOsH) + group(Cds-CdsHH) + group(Cds-CdsHH) + group(Cdd-CdsCds) + radical(Cds_P) + radical(C=C=CJ) + radical(Cds_S)"""),
)

species(
    label = '[CH2]C#CO[CH][C]=C(25027)',
    structure = SMILES('[CH2]C#CO[CH][C]=C'),
    E0 = (596.385,'kJ/mol'),
    modes = [
        HarmonicOscillator(frequencies=([1685,370,3025,407.5,1350,352.5,3000,3100,440,815,1455,1000,2100,2250,500,550,2950,3100,1380,975,1025,1650,393.486,393.573,393.583,393.592],'cm^-1')),
        HinderedRotor(inertia=(0.325852,'amu*angstrom^2'), symmetry=1, barrier=(35.7391,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.324708,'amu*angstrom^2'), symmetry=1, barrier=(35.7374,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.00108549,'amu*angstrom^2'), symmetry=1, barrier=(0.119627,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.325047,'amu*angstrom^2'), symmetry=1, barrier=(35.7411,'kJ/mol'), semiclassical=False),
    ],
    spinMultiplicity = 4,
    opticalIsomers = 1,
    molecularWeight = (93.1033,'amu'),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[1.07253,0.0594181,-5.38121e-05,2.45558e-08,-4.46326e-12,71838.1,25.6396], Tmin=(100,'K'), Tmax=(1321.65,'K')), NASAPolynomial(coeffs=[14.3699,0.0191733,-8.1364e-06,1.51594e-09,-1.05094e-13,68323.2,-42.2206], Tmin=(1321.65,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(596.385,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(266.063,'J/(mol*K)'), comment="""Thermo group additivity estimation: group(O2s-Cs(Cds-Cd)) + group(Cs-(Cds-Cds)OsHH) + group(Cs-CtHHH) + group(Cds-CdsCsH) + group(Cds-CdsHH) + group(Ct-CtCs) + group(Ct-CtOs) + radical(Cds_S) + radical(C=CCJ(O)C) + radical(Propargyl)"""),
)

species(
    label = '[CH]=C=[C]O[CH]C=C(21928)',
    structure = SMILES('[CH]=C=[C]O[CH]C=C'),
    E0 = (569.791,'kJ/mol'),
    modes = [
        HarmonicOscillator(frequencies=([3120,650,792.5,1650,3025,407.5,1350,352.5,540,610,2055,2950,3100,1380,975,1025,1650,1685,370,3010,987.5,1337.5,450,1655,331.972,331.983,332.17],'cm^-1')),
        HinderedRotor(inertia=(0.490984,'amu*angstrom^2'), symmetry=1, barrier=(38.4834,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.247295,'amu*angstrom^2'), symmetry=1, barrier=(19.3036,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(1.03731,'amu*angstrom^2'), symmetry=1, barrier=(81.1443,'kJ/mol'), semiclassical=False),
    ],
    spinMultiplicity = 4,
    opticalIsomers = 1,
    molecularWeight = (93.1033,'amu'),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[0.987093,0.0607063,-5.73201e-05,2.66754e-08,-4.46875e-12,68643.3,28.1498], Tmin=(100,'K'), Tmax=(1023.97,'K')), NASAPolynomial(coeffs=[13.9575,0.0179861,-6.38141e-06,1.0906e-09,-7.27954e-14,65570.4,-36.7664], Tmin=(1023.97,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(569.791,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(270.22,'J/(mol*K)'), comment="""Thermo group additivity estimation: group(O2s-Cs(Cds-Cd)) + group(Cs-(Cds-Cds)OsHH) + group(Cds-CdsCsH) + group(Cds-CdsOsH) + group(Cds-CdsHH) + group(Cds-CdsHH) + group(Cdd-CdsCds) + radical(C=CCJ(O)C) + radical(C=C=CJ) + radical(C=CJO)"""),
)

species(
    label = '[CH2][C]=C[O](18757)',
    structure = SMILES('[CH2][C]=C[O]'),
    E0 = (328.135,'kJ/mol'),
    modes = [
        HarmonicOscillator(frequencies=([1685,370,3000,3100,440,815,1455,1000,3010,987.5,1337.5,450,1655,180],'cm^-1')),
        HinderedRotor(inertia=(1.58249,'amu*angstrom^2'), symmetry=1, barrier=(36.3845,'kJ/mol'), semiclassical=False),
    ],
    spinMultiplicity = 4,
    opticalIsomers = 1,
    molecularWeight = (55.0553,'amu'),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[2.66329,0.0190212,1.69875e-05,-4.37428e-08,2.05267e-11,39523.3,13.9877], Tmin=(100,'K'), Tmax=(919.391,'K')), NASAPolynomial(coeffs=[12.8314,0.0017832,1.05974e-06,-2.5062e-10,1.45216e-14,36512.4,-40.4192], Tmin=(919.391,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(328.135,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(153.818,'J/(mol*K)'), comment="""Thermo group additivity estimation: group(O2s-(Cds-Cd)H) + group(Cs-(Cds-Cds)HHH) + group(Cds-CdsCsH) + group(Cds-CdsOsH) + radical(Allyl_P) + radical(C=COJ) + radical(Cds_S)"""),
)

species(
    label = '[CH]=C=[C]O[CH][C]=C(25028)',
    structure = SMILES('[CH]=C=[C]O[CH][C]=C'),
    E0 = (807.633,'kJ/mol'),
    modes = [
        HarmonicOscillator(frequencies=([1670,1700,300,440,3025,407.5,1350,352.5,540,610,2055,3120,650,792.5,1650,2950,3100,1380,975,1025,1650,333.54,333.543,333.826],'cm^-1')),
        HinderedRotor(inertia=(0.223324,'amu*angstrom^2'), symmetry=1, barrier=(17.6252,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(1.06636,'amu*angstrom^2'), symmetry=1, barrier=(84.2735,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(1.06577,'amu*angstrom^2'), symmetry=1, barrier=(84.2731,'kJ/mol'), semiclassical=False),
    ],
    spinMultiplicity = 5,
    opticalIsomers = 1,
    molecularWeight = (92.0954,'amu'),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[1.09824,0.0635674,-7.58694e-05,4.76728e-08,-1.18427e-11,97240.5,28.2278], Tmin=(100,'K'), Tmax=(986.784,'K')), NASAPolynomial(coeffs=[12.2464,0.0183767,-7.174e-06,1.26148e-09,-8.42458e-14,95040.4,-25.4069], Tmin=(986.784,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(807.633,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(245.277,'J/(mol*K)'), comment="""Thermo group additivity estimation: group(O2s-Cs(Cds-Cd)) + group(Cs-(Cds-Cds)OsHH) + group(Cds-CdsCsH) + group(Cds-CdsOsH) + group(Cds-CdsHH) + group(Cds-CdsHH) + group(Cdd-CdsCds) + radical(C=CCJ(O)C) + radical(Cds_S) + radical(C=CJO) + radical(C=C=CJ)"""),
)

species(
    label = '[CH]=[C][CH]OC=C=[CH](22637)',
    structure = SMILES('[CH]=[C][CH]OC=C=[CH]'),
    E0 = (814.984,'kJ/mol'),
    modes = [
        HarmonicOscillator(frequencies=([3115,3125,620,680,785,800,1600,1700,3025,407.5,1350,352.5,540,610,2055,3010,987.5,1337.5,450,1655,1685,370,180,180],'cm^-1')),
        HinderedRotor(inertia=(1.44189,'amu*angstrom^2'), symmetry=1, barrier=(33.1519,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(1.44444,'amu*angstrom^2'), symmetry=1, barrier=(33.2106,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(1.44257,'amu*angstrom^2'), symmetry=1, barrier=(33.1675,'kJ/mol'), semiclassical=False),
    ],
    spinMultiplicity = 5,
    opticalIsomers = 1,
    molecularWeight = (92.0954,'amu'),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[0.736887,0.0662785,-6.90866e-05,3.09145e-08,-3.59339e-12,98142.6,26.104], Tmin=(100,'K'), Tmax=(920.956,'K')), NASAPolynomial(coeffs=[16.833,0.011062,-3.08555e-06,4.61593e-10,-2.95335e-14,94554.7,-53.608], Tmin=(920.956,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(814.984,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(245.277,'J/(mol*K)'), comment="""Thermo group additivity estimation: group(O2s-Cs(Cds-Cd)) + group(Cs-(Cds-Cds)OsHH) + group(Cds-CdsCsH) + group(Cds-CdsOsH) + group(Cds-CdsHH) + group(Cds-CdsHH) + group(Cdd-CdsCds) + radical(Cds_S) + radical(C=C=CJ) + radical(C=CCJ(O)C) + radical(Cds_P)"""),
)

species(
    label = '[CH2][C]=COC1[C]=C1(25029)',
    structure = SMILES('[CH2][C]=COC1[C]=C1'),
    E0 = (749.416,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
    molecularWeight = (93.1033,'amu'),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[0.704616,0.0574163,-2.259e-05,-2.65421e-08,1.8198e-11,90266.4,24.8661], Tmin=(100,'K'), Tmax=(943.505,'K')), NASAPolynomial(coeffs=[21.5976,0.00628158,-8.2041e-07,1.35395e-10,-1.52364e-14,84657.3,-83.5472], Tmin=(943.505,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(749.416,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(270.22,'J/(mol*K)'), comment="""Thermo group additivity estimation: group(O2s-Cs(Cds-Cd)) + group(Cs-(Cds-Cds)(Cds-Cds)OsH) + group(Cs-(Cds-Cds)HHH) + group(Cds-CdsCsH) + group(Cds-CdsCsH) + group(Cds-CdsCsH) + group(Cds-CdsOsH) + ring(Cyclopropene) + radical(Cds_S) + radical(cyclopropenyl-vinyl) + radical(Allyl_P)"""),
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
    label = '[CH]C1=COC=C1[CH2](24862)',
    structure = SMILES('[CH]C1=COC=C1[CH2]'),
    E0 = (399.804,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
    molecularWeight = (93.1033,'amu'),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[1.59119,0.0352126,3.78016e-05,-7.67665e-08,3.23344e-11,48188.1,21.2291], Tmin=(100,'K'), Tmax=(964.831,'K')), NASAPolynomial(coeffs=[15.7541,0.0203078,-7.13904e-06,1.34977e-09,-1.01376e-13,43415.9,-57.1591], Tmin=(964.831,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(399.804,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(274.378,'J/(mol*K)'), comment="""Thermo group additivity estimation: group(O2s-(Cds-Cd)(Cds-Cd)) + group(Cs-(Cds-Cds)HHH) + group(Cs-(Cds-Cds)HHH) + group(Cds-Cds(Cds-Cds)Cs) + group(Cds-Cds(Cds-Cds)Cs) + group(Cds-CdsOsH) + group(Cds-CdsOsH) + ring(Furan) + radical(AllylJ2_triplet) + radical(Allyl_P)"""),
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
    label = '[CH]=[C][CH]OC#C[CH2](25030)',
    structure = SMILES('[CH]=[C][CH]OC#C[CH2]'),
    E0 = (843.482,'kJ/mol'),
    modes = [
        HarmonicOscillator(frequencies=([3120,650,792.5,1650,3025,407.5,1350,352.5,3000,3100,440,815,1455,1000,2100,2250,500,550,1685,370,438.632,439.335,440.254],'cm^-1')),
        HinderedRotor(inertia=(0.161446,'amu*angstrom^2'), symmetry=1, barrier=(22.0435,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.864417,'amu*angstrom^2'), symmetry=1, barrier=(119.627,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.867811,'amu*angstrom^2'), symmetry=1, barrier=(119.627,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.877097,'amu*angstrom^2'), symmetry=1, barrier=(119.627,'kJ/mol'), semiclassical=False),
    ],
    spinMultiplicity = 5,
    opticalIsomers = 1,
    molecularWeight = (92.0954,'amu'),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[1.1578,0.0610892,-6.41579e-05,3.40712e-08,-7.18712e-12,101551,25.7361], Tmin=(100,'K'), Tmax=(1148.48,'K')), NASAPolynomial(coeffs=[13.4706,0.0182056,-8.14901e-06,1.5595e-09,-1.10026e-13,98722.8,-35.3705], Tmin=(1148.48,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(843.482,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(241.12,'J/(mol*K)'), comment="""Thermo group additivity estimation: group(O2s-Cs(Cds-Cd)) + group(Cs-(Cds-Cds)OsHH) + group(Cs-CtHHH) + group(Cds-CdsCsH) + group(Cds-CdsHH) + group(Ct-CtCs) + group(Ct-CtOs) + radical(Propargyl) + radical(C=CCJ(O)C) + radical(Cds_S) + radical(Cds_P)"""),
)

species(
    label = '[C]=C=CO[CH][C]=C(25031)',
    structure = SMILES('[C]=C=CO[CH][C]=C'),
    E0 = (971.514,'kJ/mol'),
    modes = [
        HarmonicOscillator(frequencies=([1685,370,3025,407.5,1350,352.5,540,610,2055,3010,987.5,1337.5,450,1655,2950,3100,1380,975,1025,1650,226.748,226.795,226.836,226.867],'cm^-1')),
        HinderedRotor(inertia=(0.756551,'amu*angstrom^2'), symmetry=1, barrier=(27.6545,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.755964,'amu*angstrom^2'), symmetry=1, barrier=(27.656,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.760158,'amu*angstrom^2'), symmetry=1, barrier=(27.6588,'kJ/mol'), semiclassical=False),
    ],
    spinMultiplicity = 5,
    opticalIsomers = 1,
    molecularWeight = (92.0954,'amu'),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[0.901798,0.0664448,-7.56343e-05,4.3053e-08,-9.60899e-12,116959,25.5557], Tmin=(100,'K'), Tmax=(1096.43,'K')), NASAPolynomial(coeffs=[14.6557,0.0162678,-6.98854e-06,1.31411e-09,-9.19994e-14,113943,-42.0651], Tmin=(1096.43,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(971.514,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(245.277,'J/(mol*K)'), comment="""Thermo group additivity estimation: group(O2s-Cs(Cds-Cd)) + group(Cs-(Cds-Cds)OsHH) + group(Cds-CdsCsH) + group(Cds-CdsOsH) + group(Cds-CdsHH) + group(Cds-CdsHH) + group(Cdd-CdsCds) + radical(C=CCJ(O)C) + radical(CdCdJ2_triplet) + radical(Cds_S)"""),
)

species(
    label = '[CH]=[C]C1C[C]=CO1(25032)',
    structure = SMILES('[CH]=[C]C1C[C]=CO1'),
    E0 = (661.499,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
    molecularWeight = (93.1033,'amu'),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[1.62961,0.0328102,4.04756e-05,-9.15908e-08,4.26596e-11,79663.4,21.0333], Tmin=(100,'K'), Tmax=(903.087,'K')), NASAPolynomial(coeffs=[19.7309,0.004063,2.80327e-06,-7.22645e-10,4.84617e-14,74296.8,-76.0609], Tmin=(903.087,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(661.499,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(278.535,'J/(mol*K)'), comment="""Thermo group additivity estimation: group(O2s-Cs(Cds-Cd)) + group(Cs-(Cds-Cds)CsOsH) + group(Cs-(Cds-Cds)CsHH) + group(Cds-CdsCsH) + group(Cds-CdsCsH) + group(Cds-CdsOsH) + group(Cds-CdsHH) + ring(2,3-Dihydrofuran) + radical(Cds_S) + radical(Cds_S) + radical(Cds_P)"""),
)

species(
    label = '[CH]=C=COC#C[CH2](25033)',
    structure = SMILES('[CH]=C=COC#C[CH2]'),
    E0 = (588.422,'kJ/mol'),
    modes = [
        HarmonicOscillator(frequencies=([3120,650,792.5,1650,540,610,2055,3010,987.5,1337.5,450,1655,2100,2250,500,550,3000,3100,440,815,1455,1000,180,180],'cm^-1')),
        HinderedRotor(inertia=(1.52185,'amu*angstrom^2'), symmetry=1, barrier=(34.9904,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(1.52375,'amu*angstrom^2'), symmetry=1, barrier=(35.0341,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(1.52269,'amu*angstrom^2'), symmetry=1, barrier=(35.0096,'kJ/mol'), semiclassical=False),
    ],
    spinMultiplicity = 3,
    opticalIsomers = 1,
    molecularWeight = (92.0954,'amu'),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[0.779537,0.0624258,-6.59476e-05,3.47019e-08,-7.02934e-12,70893.9,26.019], Tmin=(100,'K'), Tmax=(1270.17,'K')), NASAPolynomial(coeffs=[16.2232,0.0116066,-3.35361e-06,4.94733e-10,-3.00758e-14,67146.9,-51.4876], Tmin=(1270.17,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(588.422,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(245.277,'J/(mol*K)'), comment="""Thermo group additivity estimation: group(O2s-(Cds-Cd)(Cds-Cd)) + group(Cs-CtHHH) + group(Cds-CdsOsH) + group(Cds-CdsHH) + group(Cdd-CdsCds) + group(Ct-CtCs) + group(Ct-CtOs) + radical(Propargyl) + radical(C=C=CJ)"""),
)

species(
    label = '[CH]=C=CO[C]=C[CH2](25034)',
    structure = SMILES('[CH]=C=CO[C]=C[CH2]'),
    E0 = (620.087,'kJ/mol'),
    modes = [
        HarmonicOscillator(frequencies=([3120,650,792.5,1650,1685,370,540,610,2055,2995,3025,975,1000,1300,1375,400,500,1630,1680,3000,3100,440,815,1455,1000,180,180],'cm^-1')),
        HinderedRotor(inertia=(1.08642,'amu*angstrom^2'), symmetry=1, barrier=(24.9789,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(2.23935,'amu*angstrom^2'), symmetry=1, barrier=(51.487,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(1.08185,'amu*angstrom^2'), symmetry=1, barrier=(24.8738,'kJ/mol'), semiclassical=False),
    ],
    spinMultiplicity = 4,
    opticalIsomers = 1,
    molecularWeight = (93.1033,'amu'),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[0.659317,0.0645191,-6.60614e-05,3.44205e-08,-6.93184e-12,74707.1,29.7314], Tmin=(100,'K'), Tmax=(1294.23,'K')), NASAPolynomial(coeffs=[15.9972,0.0140103,-3.92348e-06,5.59108e-10,-3.29286e-14,70997.1,-47.2161], Tmin=(1294.23,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(620.087,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(270.22,'J/(mol*K)'), comment="""Thermo group additivity estimation: group(O2s-(Cds-Cd)(Cds-Cd)) + group(Cs-(Cds-Cds)HHH) + group(Cds-CdsCsH) + group(Cds-CdsOsH) + group(Cds-CdsOsH) + group(Cds-CdsHH) + group(Cdd-CdsCds) + radical(Allyl_P) + radical(C=CJO) + radical(C=C=CJ)"""),
)

species(
    label = '[CH]=[C][CH]OC#CC(25035)',
    structure = SMILES('[CH]=[C][CH]OC#CC'),
    E0 = (687.075,'kJ/mol'),
    modes = [
        HarmonicOscillator(frequencies=([1685,370,2750,2800,2850,1350,1500,750,1050,1375,1000,2100,2250,500,550,3120,650,792.5,1650,3025,407.5,1350,352.5,358.231,358.236,358.238],'cm^-1')),
        HinderedRotor(inertia=(0.243009,'amu*angstrom^2'), symmetry=1, barrier=(22.1314,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.837123,'amu*angstrom^2'), symmetry=1, barrier=(76.2322,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.243031,'amu*angstrom^2'), symmetry=1, barrier=(22.1314,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.83713,'amu*angstrom^2'), symmetry=1, barrier=(76.2321,'kJ/mol'), semiclassical=False),
    ],
    spinMultiplicity = 4,
    opticalIsomers = 1,
    molecularWeight = (93.1033,'amu'),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[0.961312,0.0633128,-6.50098e-05,3.48067e-08,-7.39267e-12,82748.4,25.3827], Tmin=(100,'K'), Tmax=(1146.69,'K')), NASAPolynomial(coeffs=[13.5904,0.0192579,-7.37959e-06,1.30073e-09,-8.75688e-14,79852.1,-37.2737], Tmin=(1146.69,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(687.075,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(266.063,'J/(mol*K)'), comment="""Thermo group additivity estimation: group(O2s-Cs(Cds-Cd)) + group(Cs-(Cds-Cds)OsHH) + group(Cs-CtHHH) + group(Cds-CdsCsH) + group(Cds-CdsHH) + group(Ct-CtCs) + group(Ct-CtOs) + radical(Cds_P) + radical(Cds_S) + radical(C=CCJ(O)C)"""),
)

species(
    label = '[CH]=C=[C]OC=[C]C(25036)',
    structure = SMILES('[CH]=C=[C]OC=[C]C'),
    E0 = (706.43,'kJ/mol'),
    modes = [
        HarmonicOscillator(frequencies=([3120,650,792.5,1650,2750,2800,2850,1350,1500,750,1050,1375,1000,540,610,2055,3010,987.5,1337.5,450,1655,1670,1700,300,440,180,180],'cm^-1')),
        HinderedRotor(inertia=(0.77887,'amu*angstrom^2'), symmetry=1, barrier=(17.9077,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.778857,'amu*angstrom^2'), symmetry=1, barrier=(17.9075,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.778849,'amu*angstrom^2'), symmetry=1, barrier=(17.9073,'kJ/mol'), semiclassical=False),
    ],
    spinMultiplicity = 4,
    opticalIsomers = 1,
    molecularWeight = (93.1033,'amu'),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[0.990349,0.0673048,-8.42944e-05,5.74856e-08,-1.56039e-11,85071.2,28.739], Tmin=(100,'K'), Tmax=(903.24,'K')), NASAPolynomial(coeffs=[11.3132,0.0215908,-8.37869e-06,1.45416e-09,-9.56352e-14,83206.4,-20.0122], Tmin=(903.24,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(706.43,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(270.22,'J/(mol*K)'), comment="""Thermo group additivity estimation: group(O2s-(Cds-Cd)(Cds-Cd)) + group(Cs-(Cds-Cds)HHH) + group(Cds-CdsCsH) + group(Cds-CdsOsH) + group(Cds-CdsOsH) + group(Cds-CdsHH) + group(Cdd-CdsCds) + radical(Cds_S) + radical(C=CJO) + radical(C=C=CJ)"""),
)

species(
    label = '[CH]=C1[CH]OC=[C]C1(24777)',
    structure = SMILES('[CH]=C1[CH]OC=[C]C1'),
    E0 = (522.165,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
    molecularWeight = (93.1033,'amu'),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[1.92854,0.0222407,7.47454e-05,-1.15956e-07,4.59772e-11,62897.5,18.2446], Tmin=(100,'K'), Tmax=(971.569,'K')), NASAPolynomial(coeffs=[17.6134,0.016069,-5.89567e-06,1.25004e-09,-1.02374e-13,57093.2,-71.1594], Tmin=(971.569,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(522.165,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(282.692,'J/(mol*K)'), comment="""Thermo group additivity estimation: group(O2s-Cs(Cds-Cd)) + group(Cs-(Cds-Cds)(Cds-Cds)HH) + group(Cs-(Cds-Cds)OsHH) + group(Cds-CdsCsCs) + group(Cds-CdsCsH) + group(Cds-CdsOsH) + group(Cds-CdsHH) + ring(Cyclohexane) + radical(Cds_S) + radical(Cds_P) + radical(C=CCJ(O)C)"""),
)

species(
    label = '[CH]=[C]CO[C]=[C][CH2](25037)',
    structure = SMILES('[CH]=[C]CO[C]=[C][CH2]'),
    E0 = (1002.05,'kJ/mol'),
    modes = [
        HarmonicOscillator(frequencies=([3120,650,792.5,1650,2750,2850,1437.5,1250,1305,750,350,3000,3100,440,815,1455,1000,1670,1685,1700,300,370,440,407.275,407.276,407.276],'cm^-1')),
        HinderedRotor(inertia=(0.113542,'amu*angstrom^2'), symmetry=1, barrier=(13.3646,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.113541,'amu*angstrom^2'), symmetry=1, barrier=(13.3646,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.113541,'amu*angstrom^2'), symmetry=1, barrier=(13.3646,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.258668,'amu*angstrom^2'), symmetry=1, barrier=(30.4469,'kJ/mol'), semiclassical=False),
    ],
    spinMultiplicity = 6,
    opticalIsomers = 1,
    molecularWeight = (93.1033,'amu'),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[0.825196,0.069743,-8.26345e-05,5.0483e-08,-1.21731e-11,120633,29.7513], Tmin=(100,'K'), Tmax=(1014.99,'K')), NASAPolynomial(coeffs=[13.6368,0.0192529,-8.01723e-06,1.47248e-09,-1.01289e-13,118032,-32.2481], Tmin=(1014.99,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(1002.05,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(266.063,'J/(mol*K)'), comment="""Thermo group additivity estimation: group(O2s-Cs(Cds-Cd)) + group(Cs-(Cds-Cds)OsHH) + group(Cs-(Cds-Cds)HHH) + group(Cds-CdsCsH) + group(Cds-CdsCsH) + group(Cds-CdsOsH) + group(Cds-CdsHH) + radical(C=CJO) + radical(Cds_P) + radical(Allyl_P) + radical(Cds_S) + radical(Cds_S)"""),
)

species(
    label = '[CH]=C[CH]O[C]=[C][CH2](25038)',
    structure = SMILES('[CH]=C[CH]O[C]=[C][CH2]'),
    E0 = (875.147,'kJ/mol'),
    modes = [
        HarmonicOscillator(frequencies=([1670,1700,300,440,3025,407.5,1350,352.5,3010,987.5,1337.5,450,1655,3120,650,792.5,1650,3000,3100,440,815,1455,1000,180,377.49,861.768],'cm^-1')),
        HinderedRotor(inertia=(0.685756,'amu*angstrom^2'), symmetry=1, barrier=(15.7669,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.0298931,'amu*angstrom^2'), symmetry=1, barrier=(15.7703,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(1.55722,'amu*angstrom^2'), symmetry=1, barrier=(35.8036,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.165688,'amu*angstrom^2'), symmetry=1, barrier=(87.0849,'kJ/mol'), semiclassical=False),
    ],
    spinMultiplicity = 6,
    opticalIsomers = 1,
    molecularWeight = (93.1033,'amu'),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[1.0358,0.0632045,-6.43491e-05,3.38841e-08,-7.1252e-12,105364,29.4549], Tmin=(100,'K'), Tmax=(1150.49,'K')), NASAPolynomial(coeffs=[13.3233,0.0204834,-8.64936e-06,1.60796e-09,-1.11588e-13,102537,-31.5476], Tmin=(1150.49,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(875.147,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(266.063,'J/(mol*K)'), comment="""Thermo group additivity estimation: group(O2s-Cs(Cds-Cd)) + group(Cs-(Cds-Cds)OsHH) + group(Cs-(Cds-Cds)HHH) + group(Cds-CdsCsH) + group(Cds-CdsCsH) + group(Cds-CdsOsH) + group(Cds-CdsHH) + radical(Cds_P) + radical(Cds_S) + radical(C=CJO) + radical(Allyl_P) + radical(C=CCJ(O)C)"""),
)

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
    collisionModel = TransportData(shapeIndex=2, epsilon=(3684.58,'J/mol'), sigma=(5.9792,'angstroms'), dipoleMoment=(0,'C*m'), polarizability=(0,'angstroms^3'), rotrelaxcollnum=0, comment="""Epsilon & sigma estimated with Tc=575.52 K, Pc=39.11 bar (from Joback method)"""),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[1.29383,0.062448,-5.49929e-05,2.6954e-08,-5.66625e-12,59842.6,22.8521], Tmin=(100,'K'), Tmax=(1098.96,'K')), NASAPolynomial(coeffs=[9.25525,0.0334701,-1.54402e-05,2.96006e-09,-2.07929e-13,58092.7,-16.3084], Tmin=(1098.96,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(496.77,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(270.22,'J/(mol*K)'), comment="""Thermo group additivity estimation: group(Cs-(Cds-Cds)HHH) + group(Cd-CdCs(CO)) + group(Cds-Cds(Cds-Cds)H) + group(Cds-Cds(Cds-Cds)H) + group(Cds-O2d(Cds-Cds)H) + group(Cds-CdsHH) + radical(C=CJC=C) + radical(AllylJ2_triplet)"""),
)

species(
    label = '[CH]=[C]COC#C[CH2](24817)',
    structure = SMILES('[CH]=[C]COC#C[CH2]'),
    E0 = (732.542,'kJ/mol'),
    modes = [
        HarmonicOscillator(frequencies=([3120,650,792.5,1650,2750,2850,1437.5,1250,1305,750,350,3000,3100,440,815,1455,1000,2100,2250,500,550,1685,370,368.224,368.243,368.334],'cm^-1')),
        HinderedRotor(inertia=(0.229394,'amu*angstrom^2'), symmetry=1, barrier=(22.0698,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.229359,'amu*angstrom^2'), symmetry=1, barrier=(22.0692,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.229368,'amu*angstrom^2'), symmetry=1, barrier=(22.0691,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.229289,'amu*angstrom^2'), symmetry=1, barrier=(22.0686,'kJ/mol'), semiclassical=False),
    ],
    spinMultiplicity = 4,
    opticalIsomers = 1,
    molecularWeight = (93.1033,'amu'),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[0.732815,0.0659742,-6.80352e-05,3.48616e-08,-6.97347e-12,88227.1,26.3251], Tmin=(100,'K'), Tmax=(1223.65,'K')), NASAPolynomial(coeffs=[16.1971,0.0154226,-6.06684e-06,1.09999e-09,-7.57042e-14,84442.6,-51.4023], Tmin=(1223.65,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(732.542,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(266.063,'J/(mol*K)'), comment="""Thermo group additivity estimation: group(O2s-Cs(Cds-Cd)) + group(Cs-(Cds-Cds)OsHH) + group(Cs-CtHHH) + group(Cds-CdsCsH) + group(Cds-CdsHH) + group(Ct-CtCs) + group(Ct-CtOs) + radical(Cds_P) + radical(Propargyl) + radical(Cds_S)"""),
)

species(
    label = '[C]#CCOC=[C][CH2](25039)',
    structure = SMILES('[C]#CCOC=[C][CH2]'),
    E0 = (780.461,'kJ/mol'),
    modes = [
        HarmonicOscillator(frequencies=([1685,370,2750,2850,1437.5,1250,1305,750,350,3010,987.5,1337.5,450,1655,2175,525,3000,3100,440,815,1455,1000,263.022,263.023,263.023,263.023],'cm^-1')),
        HinderedRotor(inertia=(1.43226,'amu*angstrom^2'), symmetry=1, barrier=(70.3135,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.470501,'amu*angstrom^2'), symmetry=1, barrier=(23.098,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.470504,'amu*angstrom^2'), symmetry=1, barrier=(23.0981,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(1.43226,'amu*angstrom^2'), symmetry=1, barrier=(70.3134,'kJ/mol'), semiclassical=False),
    ],
    spinMultiplicity = 4,
    opticalIsomers = 1,
    molecularWeight = (93.1033,'amu'),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[0.0708013,0.0740992,-8.39524e-05,4.63933e-08,-9.57904e-12,94019.9,26.8011], Tmin=(100,'K'), Tmax=(1359.83,'K')), NASAPolynomial(coeffs=[18.6935,0.0088777,-4.89448e-07,-1.72073e-10,2.00174e-14,89920.6,-65.2166], Tmin=(1359.83,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(780.461,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(266.063,'J/(mol*K)'), comment="""Thermo group additivity estimation: group(O2s-Cs(Cds-Cd)) + group(Cs-CtOsHH) + group(Cs-(Cds-Cds)HHH) + group(Cds-CdsCsH) + group(Cds-CdsOsH) + group(Ct-CtCs) + group(Ct-CtH) + radical(Acetyl) + radical(Cds_S) + radical(Allyl_P)"""),
)

species(
    label = '[C]=C=CO[CH]C=C(25040)',
    structure = SMILES('[C]=C=CO[CH]C=C'),
    E0 = (733.672,'kJ/mol'),
    modes = [
        HarmonicOscillator(frequencies=([3025,407.5,1350,352.5,540,610,2055,2995,3025,975,1000,1300,1375,400,500,1630,1680,2950,3100,1380,975,1025,1650,199.864,199.864,199.864,199.864],'cm^-1')),
        HinderedRotor(inertia=(1.02761,'amu*angstrom^2'), symmetry=1, barrier=(29.1288,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(1.02761,'amu*angstrom^2'), symmetry=1, barrier=(29.1288,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(1.02761,'amu*angstrom^2'), symmetry=1, barrier=(29.1288,'kJ/mol'), semiclassical=False),
    ],
    spinMultiplicity = 4,
    opticalIsomers = 1,
    molecularWeight = (93.1033,'amu'),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[0.581421,0.065945,-6.48366e-05,3.13752e-08,-5.89451e-12,88371.2,26.2353], Tmin=(100,'K'), Tmax=(1306,'K')), NASAPolynomial(coeffs=[17.5523,0.0139665,-5.13667e-06,9.00274e-10,-6.08462e-14,83938.5,-60.17], Tmin=(1306,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(733.672,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(270.22,'J/(mol*K)'), comment="""Thermo group additivity estimation: group(O2s-Cs(Cds-Cd)) + group(Cs-(Cds-Cds)OsHH) + group(Cds-CdsCsH) + group(Cds-CdsOsH) + group(Cds-CdsHH) + group(Cds-CdsHH) + group(Cdd-CdsCds) + radical(CdCdJ2_triplet) + radical(C=CCJ(O)C)"""),
)

species(
    label = '[C]#C[CH]OC=[C]C(25041)',
    structure = SMILES('[C]#C[CH]OC=[C]C'),
    E0 = (822.889,'kJ/mol'),
    modes = [
        HarmonicOscillator(frequencies=([1685,370,2750,2800,2850,1350,1500,750,1050,1375,1000,3010,987.5,1337.5,450,1655,2175,525,3025,407.5,1350,352.5,180,180,180,180],'cm^-1')),
        HinderedRotor(inertia=(0.917342,'amu*angstrom^2'), symmetry=1, barrier=(21.0915,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.917158,'amu*angstrom^2'), symmetry=1, barrier=(21.0873,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.916531,'amu*angstrom^2'), symmetry=1, barrier=(21.0728,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.916837,'amu*angstrom^2'), symmetry=1, barrier=(21.0799,'kJ/mol'), semiclassical=False),
    ],
    spinMultiplicity = 4,
    opticalIsomers = 1,
    molecularWeight = (93.1033,'amu'),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[0.287756,0.0788536,-9.31187e-05,4.61292e-08,-5.40051e-12,99107.3,25.2322], Tmin=(100,'K'), Tmax=(823.361,'K')), NASAPolynomial(coeffs=[18.5977,0.00911873,-1.08609e-06,-4.10131e-11,1.10687e-14,95440.8,-63.4999], Tmin=(823.361,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(822.889,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(266.063,'J/(mol*K)'), comment="""Thermo group additivity estimation: group(O2s-Cs(Cds-Cd)) + group(Cs-CtOsHH) + group(Cs-(Cds-Cds)HHH) + group(Cds-CdsCsH) + group(Cds-CdsOsH) + group(Ct-CtCs) + group(Ct-CtH) + radical(Cds_S) + radical(Acetyl) + radical(CCsJOC(O))"""),
)

species(
    label = '[CH2]C1[CH]OC=[C]C=1(24893)',
    structure = SMILES('[CH2]C1[CH]OC=[C]C=1'),
    E0 = (337.66,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
    molecularWeight = (93.1033,'amu'),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[1.76851,0.028721,5.53325e-05,-1.01712e-07,4.39632e-11,40710.3,17.5653], Tmin=(100,'K'), Tmax=(925.905,'K')), NASAPolynomial(coeffs=[17.8755,0.0114432,-1.4142e-06,1.58528e-10,-1.62543e-14,35485.5,-71.0092], Tmin=(925.905,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(337.66,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(278.535,'J/(mol*K)'), comment="""Thermo group additivity estimation: group(O2s-Cs(Cds-Cd)) + group(Cs-(Cds-Cds)OsHH) + group(Cs-(Cds-Cds)HHH) + group(Cds-CdsCsCs) + group(Cds-Cds(Cds-Cds)H) + group(Cds-Cds(Cds-Cds)H) + group(Cds-CdsOsH) + ring(1,3-Cyclohexadiene) + radical(C=CJC=C) + radical(C=CC=CCJ) + radical(C=CCJ(O)C)"""),
)

species(
    label = '[C]1[CH]OC=[C]CC=1(24811)',
    structure = SMILES('[C]1[CH]OC=[C]CC=1'),
    E0 = (546.109,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
    molecularWeight = (93.1033,'amu'),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[2.04399,-0.0153176,0.000246569,-3.57549e-07,1.50985e-10,65807.2,22.6827], Tmin=(100,'K'), Tmax=(899.037,'K')), NASAPolynomial(coeffs=[41.4069,-0.0372165,2.74413e-05,-5.47461e-09,3.63508e-13,52536.7,-197.472], Tmin=(899.037,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(546.109,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(282.692,'J/(mol*K)'), comment="""Thermo group additivity estimation: group(O2s-Cs(Cds-Cd)) + group(Cs-(Cds-Cds)(Cds-Cds)HH) + group(Cs-(Cds-Cds)OsHH) + group(Cds-CdsCsH) + group(Cds-CdsCsH) + group(Cds-CdsCsH) + group(Cds-CdsOsH) + ring(Cycloheptane) + radical(Cds_S) + radical(Cds_S) + radical(C=CCJ(O)C)"""),
)

species(
    label = 'C#CC1OC=C1[CH2](24933)',
    structure = SMILES('C#CC1OC=C1[CH2]'),
    E0 = (341.142,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
    molecularWeight = (93.1033,'amu'),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[1.06135,0.0390468,4.71296e-05,-1.11604e-07,5.21815e-11,41159.6,19.5359], Tmin=(100,'K'), Tmax=(914.45,'K')), NASAPolynomial(coeffs=[25.9271,-0.0026836,5.61716e-06,-1.17209e-09,7.37572e-14,33809,-113.529], Tmin=(914.45,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(341.142,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(274.378,'J/(mol*K)'), comment="""Thermo group additivity estimation: group(O2s-Cs(Cds-Cd)) + group(Cs-(Cds-Cds)(Cds-Cds)OsH) + group(Cs-(Cds-Cds)HHH) + group(Cds-CdsCsCs) + group(Cds-CdsOsH) + group(Ct-CtCs) + group(Ct-CtH) + ring(Cyclobutene) + radical(Allyl_P)"""),
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
    label = '[CH]OC=[C][CH2](19096)',
    structure = SMILES('[CH]OC=[C][CH2]'),
    E0 = (666.345,'kJ/mol'),
    modes = [
        HarmonicOscillator(frequencies=([1685,370,3010,987.5,1337.5,450,1655,3000,3100,440,815,1455,1000,405.585,405.586,405.587,405.596,405.707],'cm^-1')),
        HinderedRotor(inertia=(0.00102476,'amu*angstrom^2'), symmetry=1, barrier=(0.119627,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.169574,'amu*angstrom^2'), symmetry=1, barrier=(19.7978,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.169596,'amu*angstrom^2'), symmetry=1, barrier=(19.7979,'kJ/mol'), semiclassical=False),
    ],
    spinMultiplicity = 5,
    opticalIsomers = 1,
    molecularWeight = (68.074,'amu'),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[1.50517,0.0438729,-2.19563e-05,-1.72816e-08,1.38466e-11,80242.8,18.1997], Tmin=(100,'K'), Tmax=(920.957,'K')), NASAPolynomial(coeffs=[18.1259,0.000529453,1.65684e-06,-3.65253e-10,2.24492e-14,75958.1,-67.2588], Tmin=(920.957,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(666.345,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(195.39,'J/(mol*K)'), comment="""Thermo group additivity estimation: group(O2s-Cs(Cds-Cd)) + group(Cs-(Cds-Cds)HHH) + group(Cs-OsHHH) + group(Cds-CdsCsH) + group(Cds-CdsOsH) + radical(Allyl_P) + radical(CH2_triplet) + radical(Cds_S)"""),
)

species(
    label = '[C]=C=COC[C]=C(24330)',
    structure = SMILES('[C]=C=COC[C]=C'),
    E0 = (860.574,'kJ/mol'),
    modes = [
        HarmonicOscillator(frequencies=([1685,370,2750,2850,1437.5,1250,1305,750,350,540,610,2055,3010,987.5,1337.5,450,1655,2950,3100,1380,975,1025,1650,180,180,180,180],'cm^-1')),
        HinderedRotor(inertia=(0.973582,'amu*angstrom^2'), symmetry=1, barrier=(22.3846,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.974701,'amu*angstrom^2'), symmetry=1, barrier=(22.4103,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.974363,'amu*angstrom^2'), symmetry=1, barrier=(22.4025,'kJ/mol'), semiclassical=False),
    ],
    spinMultiplicity = 4,
    opticalIsomers = 1,
    molecularWeight = (93.1033,'amu'),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[0.499864,0.0710787,-7.87272e-05,4.29502e-08,-9.06363e-12,103634,26.0606], Tmin=(100,'K'), Tmax=(1168.17,'K')), NASAPolynomial(coeffs=[17.2677,0.0136633,-5.0026e-06,8.76177e-10,-5.93974e-14,99716.7,-57.4407], Tmin=(1168.17,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(860.574,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(270.22,'J/(mol*K)'), comment="""Thermo group additivity estimation: group(O2s-Cs(Cds-Cd)) + group(Cs-(Cds-Cds)OsHH) + group(Cds-CdsCsH) + group(Cds-CdsOsH) + group(Cds-CdsHH) + group(Cds-CdsHH) + group(Cdd-CdsCds) + radical(CdCdJ2_triplet) + radical(Cds_S)"""),
)

species(
    label = 'C=[C]C1[CH][C]=CO1(24991)',
    structure = SMILES('C=[C]C1[CH][C]=CO1'),
    E0 = (484.546,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
    molecularWeight = (93.1033,'amu'),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[1.94729,0.0152531,0.00010658,-1.68565e-07,7.17185e-11,58379.2,20.2879], Tmin=(100,'K'), Tmax=(911.171,'K')), NASAPolynomial(coeffs=[23.5777,-0.00177502,6.32414e-06,-1.34809e-09,8.53346e-14,51202.5,-99.8058], Tmin=(911.171,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(484.546,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(278.535,'J/(mol*K)'), comment="""Thermo group additivity estimation: group(O2s-Cs(Cds-Cd)) + group(Cs-(Cds-Cds)CsOsH) + group(Cs-(Cds-Cds)CsHH) + group(Cds-CdsCsH) + group(Cds-CdsCsH) + group(Cds-CdsOsH) + group(Cds-CdsHH) + ring(2,3-Dihydrofuran) + radical(Cds_S) + radical(C=CCJC(O)C=C) + radical(Cds_S)"""),
)

species(
    label = 'C#CC1OC1[C]=C(25042)',
    structure = SMILES('C#CC1OC1[C]=C'),
    E0 = (471.27,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
    molecularWeight = (93.1033,'amu'),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[-0.211414,0.062905,-5.9304e-05,2.83306e-08,-4.91502e-12,56858.2,28.453], Tmin=(100,'K'), Tmax=(1771.75,'K')), NASAPolynomial(coeffs=[13.308,0.0124625,2.66662e-07,-4.30381e-10,3.86824e-14,55194.2,-35.6795], Tmin=(1771.75,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(471.27,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(274.378,'J/(mol*K)'), comment="""Thermo group additivity estimation: group(O2s-CsCs) + group(Cs-(Cds-Cds)CsOsH) + group(Cs-(Cds-Cds)CsOsH) + group(Cds-CdsCsH) + group(Cds-CdsHH) + group(Ct-CtCs) + group(Ct-CtH) + ring(Ethylene_oxide) + radical(Cds_S)"""),
)

species(
    label = '[CH]=C[CH]OC#C[CH2](25043)',
    structure = SMILES('[CH]=C[CH]OC#C[CH2]'),
    E0 = (605.64,'kJ/mol'),
    modes = [
        HarmonicOscillator(frequencies=([3120,650,792.5,1650,3025,407.5,1350,352.5,3010,987.5,1337.5,450,1655,2100,2250,500,550,3000,3100,440,815,1455,1000,461.541,461.542,461.542],'cm^-1')),
        HinderedRotor(inertia=(0.182847,'amu*angstrom^2'), symmetry=1, barrier=(27.6399,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.684535,'amu*angstrom^2'), symmetry=1, barrier=(103.477,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.182847,'amu*angstrom^2'), symmetry=1, barrier=(27.6399,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.182847,'amu*angstrom^2'), symmetry=1, barrier=(27.6399,'kJ/mol'), semiclassical=False),
    ],
    spinMultiplicity = 4,
    opticalIsomers = 1,
    molecularWeight = (93.1033,'amu'),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[0.910496,0.0597665,-5.06814e-05,1.92334e-08,-2.26784e-12,72959.8,26.1512], Tmin=(100,'K'), Tmax=(1113.82,'K')), NASAPolynomial(coeffs=[15.6485,0.0170624,-6.93966e-06,1.2932e-09,-9.08497e-14,69042.5,-49.3869], Tmin=(1113.82,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(605.64,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(266.063,'J/(mol*K)'), comment="""Thermo group additivity estimation: group(O2s-Cs(Cds-Cd)) + group(Cs-(Cds-Cds)OsHH) + group(Cs-CtHHH) + group(Cds-CdsCsH) + group(Cds-CdsHH) + group(Ct-CtCs) + group(Ct-CtOs) + radical(Cds_P) + radical(C=CCJ(O)C) + radical(Propargyl)"""),
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
    label = '[CH]=[C][CH]O[C]=C[CH2](21941)',
    structure = SMILES('[CH]=[C][CH]O[C]=C[CH2]'),
    E0 = (875.147,'kJ/mol'),
    modes = [
        HarmonicOscillator(frequencies=([1670,1700,300,440,3025,407.5,1350,352.5,3010,987.5,1337.5,450,1655,3120,650,792.5,1650,3000,3100,440,815,1455,1000,180,180,861.58],'cm^-1')),
        HinderedRotor(inertia=(0.0299371,'amu*angstrom^2'), symmetry=1, barrier=(15.7686,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.155953,'amu*angstrom^2'), symmetry=1, barrier=(15.7689,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(1.55724,'amu*angstrom^2'), symmetry=1, barrier=(35.804,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.165361,'amu*angstrom^2'), symmetry=1, barrier=(87.091,'kJ/mol'), semiclassical=False),
    ],
    spinMultiplicity = 6,
    opticalIsomers = 1,
    molecularWeight = (93.1033,'amu'),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[1.0358,0.0632045,-6.43492e-05,3.38842e-08,-7.12521e-12,105364,29.4549], Tmin=(100,'K'), Tmax=(1150.48,'K')), NASAPolynomial(coeffs=[13.3233,0.0204834,-8.64937e-06,1.60796e-09,-1.11588e-13,102537,-31.5476], Tmin=(1150.48,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(875.147,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(266.063,'J/(mol*K)'), comment="""Thermo group additivity estimation: group(O2s-Cs(Cds-Cd)) + group(Cs-(Cds-Cds)OsHH) + group(Cs-(Cds-Cds)HHH) + group(Cds-CdsCsH) + group(Cds-CdsCsH) + group(Cds-CdsOsH) + group(Cds-CdsHH) + radical(C=CCJ(O)C) + radical(Allyl_P) + radical(C=CJO) + radical(Cds_S) + radical(Cds_P)"""),
)

species(
    label = '[CH]C=CO[CH][C]=[CH](21936)',
    structure = SMILES('[CH]C=CO[CH][C]=[CH]'),
    E0 = (854.588,'kJ/mol'),
    modes = [
        HarmonicOscillator(frequencies=([3120,650,792.5,1650,3025,407.5,1350,352.5,2995,3025,975,1000,1300,1375,400,500,1630,1680,1685,370,200,800,1000,1200,1400,1600],'cm^-1')),
        HinderedRotor(inertia=(0.156089,'amu*angstrom^2'), symmetry=1, barrier=(3.5888,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.156089,'amu*angstrom^2'), symmetry=1, barrier=(3.5888,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.156089,'amu*angstrom^2'), symmetry=1, barrier=(3.5888,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.156089,'amu*angstrom^2'), symmetry=1, barrier=(3.5888,'kJ/mol'), semiclassical=False),
    ],
    spinMultiplicity = 6,
    opticalIsomers = 1,
    molecularWeight = (93.1033,'amu'),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[0.615174,0.0663702,-5.44861e-05,1.95288e-08,-1.71177e-12,102912,28.6628], Tmin=(100,'K'), Tmax=(1048.69,'K')), NASAPolynomial(coeffs=[15.5714,0.0217291,-8.37861e-06,1.49852e-09,-1.0277e-13,99092.9,-47.456], Tmin=(1048.69,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(854.588,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(266.063,'J/(mol*K)'), comment="""Thermo group additivity estimation: group(O2s-Cs(Cds-Cd)) + group(Cs-(Cds-Cds)OsHH) + group(Cs-(Cds-Cds)HHH) + group(Cds-CdsCsH) + group(Cds-CdsCsH) + group(Cds-CdsOsH) + group(Cds-CdsHH) + radical(C=CCJ(O)C) + radical(Cds_P) + radical(AllylJ2_triplet) + radical(Cds_S)"""),
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
    label = '[CH]=C=[C]OC=C=C(25044)',
    structure = SMILES('[CH]=C=[C]OC=C=C'),
    E0 = (645.192,'kJ/mol'),
    modes = [
        HarmonicOscillator(frequencies=([1685,370,2950,3100,1380,975,1025,1650,3010,987.5,1337.5,450,1655,3120,650,792.5,1650,540,563.333,586.667,610,1970,2140,180,180],'cm^-1')),
        HinderedRotor(inertia=(1.35741,'amu*angstrom^2'), symmetry=1, barrier=(31.2095,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(1.35722,'amu*angstrom^2'), symmetry=1, barrier=(31.2053,'kJ/mol'), semiclassical=False),
    ],
    spinMultiplicity = 3,
    opticalIsomers = 1,
    molecularWeight = (92.0954,'amu'),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[0.987529,0.0632695,-7.19116e-05,4.11471e-08,-8.92407e-12,77709.8,27.3262], Tmin=(100,'K'), Tmax=(950.295,'K')), NASAPolynomial(coeffs=[13.8998,0.0148926,-4.97987e-06,8.06922e-10,-5.17363e-14,74986,-35.7288], Tmin=(950.295,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(645.192,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(249.434,'J/(mol*K)'), comment="""Thermo group additivity estimation: group(O2s-(Cds-Cd)(Cds-Cd)) + group(Cds-CdsOsH) + group(Cds-CdsOsH) + group(Cds-CdsHH) + group(Cds-CdsHH) + group(Cdd-CdsCds) + group(Cdd-CdsCds) + radical(C=CJO) + radical(C=C=CJ)"""),
)

species(
    label = '[CH]C=[C]OC=C=C(25045)',
    structure = SMILES('[CH]C=[C]OC=C=C'),
    E0 = (684.796,'kJ/mol'),
    modes = [
        HarmonicOscillator(frequencies=([1685,370,2950,3100,1380,975,1025,1650,2995,3025,975,1000,1300,1375,400,500,1630,1680,540,610,2055,180,180,180,180,180,180],'cm^-1')),
        HinderedRotor(inertia=(2.15901,'amu*angstrom^2'), symmetry=1, barrier=(49.6399,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(2.15767,'amu*angstrom^2'), symmetry=1, barrier=(49.609,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(2.15706,'amu*angstrom^2'), symmetry=1, barrier=(49.595,'kJ/mol'), semiclassical=False),
    ],
    spinMultiplicity = 4,
    opticalIsomers = 1,
    molecularWeight = (93.1033,'amu'),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[0.956627,0.0622604,-5.32716e-05,2.41602e-08,-4.46082e-12,82475.3,28.8675], Tmin=(100,'K'), Tmax=(1288.58,'K')), NASAPolynomial(coeffs=[12.9229,0.0251149,-1.00316e-05,1.78935e-09,-1.20604e-13,79391.4,-31.8967], Tmin=(1288.58,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(684.796,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(270.22,'J/(mol*K)'), comment="""Thermo group additivity estimation: group(O2s-(Cds-Cd)(Cds-Cd)) + group(Cs-(Cds-Cds)HHH) + group(Cds-CdsCsH) + group(Cds-CdsOsH) + group(Cds-CdsOsH) + group(Cds-CdsHH) + group(Cdd-CdsCds) + radical(C=CJO) + radical(AllylJ2_triplet)"""),
)

species(
    label = '[CH][C]=[C]OC[C]=C(25046)',
    structure = SMILES('[CH][C]=[C]OC[C]=C'),
    E0 = (974.139,'kJ/mol'),
    modes = [
        HarmonicOscillator(frequencies=([1670,1685,1700,300,370,440,2750,2850,1437.5,1250,1305,750,350,2950,3100,1380,975,1025,1650,200,800,960,1120,1280,1440,1600],'cm^-1')),
        HinderedRotor(inertia=(0.156089,'amu*angstrom^2'), symmetry=1, barrier=(3.5888,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.156089,'amu*angstrom^2'), symmetry=1, barrier=(3.5888,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.156089,'amu*angstrom^2'), symmetry=1, barrier=(3.5888,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.156089,'amu*angstrom^2'), symmetry=1, barrier=(3.5888,'kJ/mol'), semiclassical=False),
    ],
    spinMultiplicity = 6,
    opticalIsomers = 1,
    molecularWeight = (93.1033,'amu'),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[0.986093,0.0675792,-7.0281e-05,4.06626e-08,-9.69209e-12,117269,29.6014], Tmin=(100,'K'), Tmax=(1006.24,'K')), NASAPolynomial(coeffs=[10.7611,0.0287215,-1.2356e-05,2.2855e-09,-1.57315e-13,115302,-17.6184], Tmin=(1006.24,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(974.139,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(266.063,'J/(mol*K)'), comment="""Thermo group additivity estimation: group(O2s-Cs(Cds-Cd)) + group(Cs-(Cds-Cds)OsHH) + group(Cs-(Cds-Cds)HHH) + group(Cds-CdsCsH) + group(Cds-CdsCsH) + group(Cds-CdsOsH) + group(Cds-CdsHH) + radical(AllylJ2_triplet) + radical(C=CJO) + radical(Cds_S) + radical(Cds_S)"""),
)

species(
    label = '[CH][C]=[C]O[CH]C=C(21935)',
    structure = SMILES('[CH][C]=[C]O[CH]C=C'),
    E0 = (847.236,'kJ/mol'),
    modes = [
        HarmonicOscillator(frequencies=([1670,1700,300,440,3025,407.5,1350,352.5,3010,987.5,1337.5,450,1655,2950,3100,1380,975,1025,1650,200,800,960,1120,1280,1440,1600],'cm^-1')),
        HinderedRotor(inertia=(0.156089,'amu*angstrom^2'), symmetry=1, barrier=(3.5888,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.156089,'amu*angstrom^2'), symmetry=1, barrier=(3.5888,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.156089,'amu*angstrom^2'), symmetry=1, barrier=(3.5888,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.156089,'amu*angstrom^2'), symmetry=1, barrier=(3.5888,'kJ/mol'), semiclassical=False),
    ],
    spinMultiplicity = 6,
    opticalIsomers = 1,
    molecularWeight = (93.1033,'amu'),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[1.21044,0.0609227,-5.17518e-05,2.39185e-08,-4.6308e-12,102000,29.2523], Tmin=(100,'K'), Tmax=(1205.86,'K')), NASAPolynomial(coeffs=[10.6562,0.0295895,-1.27756e-05,2.37011e-09,-1.63352e-13,99721.6,-18.0863], Tmin=(1205.86,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(847.236,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(266.063,'J/(mol*K)'), comment="""Thermo group additivity estimation: group(O2s-Cs(Cds-Cd)) + group(Cs-(Cds-Cds)OsHH) + group(Cs-(Cds-Cds)HHH) + group(Cds-CdsCsH) + group(Cds-CdsCsH) + group(Cds-CdsOsH) + group(Cds-CdsHH) + radical(C=CJO) + radical(Cds_S) + radical(AllylJ2_triplet) + radical(C=CCJ(O)C)"""),
)

species(
    label = '[CH]=C=[C]O[CH][C]C(25047)',
    structure = SMILES('[CH]=C=[C]O[CH][C]C'),
    E0 = (967.655,'kJ/mol'),
    modes = [
        HarmonicOscillator(frequencies=([1685,370,2750,2800,2850,1350,1500,750,1050,1375,1000,540,610,2055,3120,650,792.5,1650,3025,407.5,1350,352.5,180,180,180,180],'cm^-1')),
        HinderedRotor(inertia=(0.982,'amu*angstrom^2'), symmetry=1, barrier=(22.5781,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.984666,'amu*angstrom^2'), symmetry=1, barrier=(22.6394,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.983875,'amu*angstrom^2'), symmetry=1, barrier=(22.6212,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.981741,'amu*angstrom^2'), symmetry=1, barrier=(22.5722,'kJ/mol'), semiclassical=False),
    ],
    spinMultiplicity = 6,
    opticalIsomers = 1,
    molecularWeight = (93.1033,'amu'),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[-0.135759,0.085683,-0.000112087,7.04317e-08,-1.67669e-11,116536,29.7625], Tmin=(100,'K'), Tmax=(1095.6,'K')), NASAPolynomial(coeffs=[19.8133,0.00901507,-1.87064e-06,1.70905e-10,-5.48368e-15,112394,-67.2516], Tmin=(1095.6,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(967.655,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(266.063,'J/(mol*K)'), comment="""Thermo group additivity estimation: group(O2s-Cs(Cds-Cd)) + group(Cs-CsCsHH) + group(Cs-CsOsHH) + group(Cs-CsHHH) + group(Cds-CdsOsH) + group(Cds-CdsHH) + group(Cdd-CdsCds) + radical(C=C=CJ) + radical(CCsJOC(O)) + radical(CCJ2_triplet) + radical(C=CJO)"""),
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
    E0 = (567.888,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS2',
    E0 = (790.257,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS3',
    E0 = (682.194,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS4',
    E0 = (786.839,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS5',
    E0 = (636.9,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS6',
    E0 = (732.011,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS7',
    E0 = (838.67,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS8',
    E0 = (849.23,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS9',
    E0 = (762.506,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS10',
    E0 = (782.352,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS11',
    E0 = (819.816,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS12',
    E0 = (1019.44,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS13',
    E0 = (1026.79,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS14',
    E0 = (749.416,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS15',
    E0 = (694.527,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS16',
    E0 = (611.205,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS17',
    E0 = (889.691,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS18',
    E0 = (1055.29,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS19',
    E0 = (1183.32,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS20',
    E0 = (679.315,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS21',
    E0 = (815.336,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS22',
    E0 = (644.942,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS23',
    E0 = (824.224,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS24',
    E0 = (848.996,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS25',
    E0 = (739.47,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS26',
    E0 = (647.892,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS27',
    E0 = (1065.45,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS28',
    E0 = (900.12,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS29',
    E0 = (881.688,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS30',
    E0 = (721.538,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS31',
    E0 = (891.952,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS32',
    E0 = (927.738,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS33',
    E0 = (808.116,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS34',
    E0 = (906.465,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS35',
    E0 = (598.921,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS36',
    E0 = (675.112,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS37',
    E0 = (575.42,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS38',
    E0 = (1252.6,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS39',
    E0 = (912.706,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS40',
    E0 = (623.696,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS41',
    E0 = (575.749,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS42',
    E0 = (755.712,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS43',
    E0 = (875.801,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS44',
    E0 = (897.448,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS45',
    E0 = (877.449,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS46',
    E0 = (982.311,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS47',
    E0 = (872.107,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS48',
    E0 = (888.933,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS49',
    E0 = (1052.39,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS50',
    E0 = (872.209,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS51',
    E0 = (992.628,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS52',
    E0 = (881.688,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

reaction(
    label = 'reaction1',
    reactants = ['[CH]=C=CO[CH][C]=C(21927)'],
    products = ['C#CC=O(21959)', 'C#C[CH2](17441)'],
    transitionState = 'TS1',
    kinetics = Arrhenius(A=(5e+12,'s^-1'), n=0, Ea=(0,'kJ/mol'), T0=(1,'K'), Tmin=(300,'K'), Tmax=(1500,'K'), comment="""Exact match found for rate rule [RJJ]
Euclidian distance = 0
family: 1,4_Linear_birad_scission"""),
)

reaction(
    label = 'reaction2',
    reactants = ['[CH]=C=CO[CH][C]=C(21927)'],
    products = ['[CH]=[C]C1OC1[C]=C(25025)'],
    transitionState = 'TS2',
    kinetics = Arrhenius(A=(7.18499e+11,'s^-1'), n=-0.0609598, Ea=(222.369,'kJ/mol'), T0=(1,'K'), comment="""Estimated using an average for rate rule [R4_S_D;doublebond_intra;radadd_intra_cs]
Euclidian distance = 0
family: Intra_R_Add_Exocyclic
Ea raised from 220.6 to 222.4 kJ/mol to match endothermicity of reaction."""),
)

reaction(
    label = 'reaction3',
    reactants = ['[CH]=C=CO[CH][C]=C(21927)'],
    products = ['[CH]=[C]C1OC=C1[CH2](25026)'],
    transitionState = 'TS3',
    kinetics = Arrhenius(A=(9.65748e+09,'s^-1'), n=0.409794, Ea=(114.306,'kJ/mol'), T0=(1,'K'), comment="""Estimated using template [R5;doublebond_intra;radadd_intra] for rate rule [R5;doublebond_intra;radadd_intra_cddouble]
Euclidian distance = 1.0
family: Intra_R_Add_Exocyclic"""),
)

reaction(
    label = 'reaction4',
    reactants = ['H(8)', '[CH]=C=COC=C=[CH](22632)'],
    products = ['[CH]=C=CO[CH][C]=C(21927)'],
    transitionState = 'TS4',
    kinetics = Arrhenius(A=(8556.54,'m^3/(mol*s)'), n=1.383, Ea=(15.1097,'kJ/mol'), T0=(1,'K'), comment="""Estimated using an average for rate rule [Ct_Ct;HJ]
Euclidian distance = 0
Multiplied by reaction path degeneracy 2.0
family: R_Addition_MultipleBond"""),
)

reaction(
    label = 'reaction5',
    reactants = ['C=C=C[O](12570)', '[CH]=C=[CH](18734)'],
    products = ['[CH]=C=CO[CH][C]=C(21927)'],
    transitionState = 'TS5',
    kinetics = Arrhenius(A=(3.63349,'m^3/(mol*s)'), n=2.00263, Ea=(29.8204,'kJ/mol'), T0=(1,'K'), comment="""Estimated using average of templates [R_R;CJ] + [Od_R;YJ] for rate rule [Od_R;CJ]
Euclidian distance = 1.0
Multiplied by reaction path degeneracy 2.0
family: R_Addition_MultipleBond"""),
)

reaction(
    label = 'reaction6',
    reactants = ['[CH]=C=CO[CH]C=[CH](24816)'],
    products = ['[CH]=C=CO[CH][C]=C(21927)'],
    transitionState = 'TS6',
    kinetics = Arrhenius(A=(3.3587e+07,'s^-1'), n=1.74167, Ea=(154.868,'kJ/mol'), T0=(1,'K'), comment="""Estimated using an average for rate rule [R2H_D;Cd_rad_out_singleH;XH_out]
Euclidian distance = 0
family: intra_H_migration"""),
)

reaction(
    label = 'reaction7',
    reactants = ['[CH]=C=[C]OC[C]=C(24326)'],
    products = ['[CH]=C=CO[CH][C]=C(21927)'],
    transitionState = 'TS7',
    kinetics = Arrhenius(A=(4.823e+09,'s^-1'), n=1.00333, Ea=(141.977,'kJ/mol'), T0=(1,'K'), comment="""Estimated using template [R3H_SS;Cd_rad_out_double;XH_out] for rate rule [R3H_SS_O;Cd_rad_out_double;XH_out]
Euclidian distance = 1.0
Multiplied by reaction path degeneracy 2.0
family: intra_H_migration"""),
)

reaction(
    label = 'reaction8',
    reactants = ['[CH]=[C]COC=C=[CH](22625)'],
    products = ['[CH]=C=CO[CH][C]=C(21927)'],
    transitionState = 'TS8',
    kinetics = Arrhenius(A=(1.846e+10,'s^-1'), n=0.74, Ea=(145.185,'kJ/mol'), T0=(1,'K'), comment="""Estimated using template [R3Hall;Cd_rad_out_singleH;Cs_H_out_1H] for rate rule [R3HJ;Cd_rad_out_singleH;Cs_H_out_H/NonDeO]
Euclidian distance = 1.41421356237
Multiplied by reaction path degeneracy 2.0
family: intra_H_migration"""),
)

reaction(
    label = 'reaction9',
    reactants = ['[CH2]C#CO[CH][C]=C(25027)'],
    products = ['[CH]=C=CO[CH][C]=C(21927)'],
    transitionState = 'TS9',
    kinetics = Arrhenius(A=(1.65028e+09,'s^-1'), n=1.32317, Ea=(166.12,'kJ/mol'), T0=(1,'K'), comment="""Estimated using average of templates [RnH;Cd_rad_out_singleNd;Cd_H_out_singleH] + [R3H;Cd_rad_out_singleNd;XH_out] for rate rule [R3H;Cd_rad_out_singleNd;Cd_H_out_singleH]
Euclidian distance = 2.0
Multiplied by reaction path degeneracy 2.0
family: intra_H_migration"""),
)

reaction(
    label = 'reaction10',
    reactants = ['[CH]=C=CO[CH][C]=C(21927)'],
    products = ['[CH]=C=[C]O[CH]C=C(21928)'],
    transitionState = 'TS10',
    kinetics = Arrhenius(A=(1.806e+09,'s^-1'), n=1.172, Ea=(214.463,'kJ/mol'), T0=(1,'K'), comment="""Estimated using template [RnH;Cd_rad_out_Cd;Cd_H_out_doubleC] for rate rule [R4HJ_1;Cd_rad_out_Cd;Cd_H_out_doubleC]
Euclidian distance = 2.0
family: intra_H_migration"""),
)

reaction(
    label = 'reaction11',
    reactants = ['[CH2][C]=C[O](18757)', '[CH]=C=[CH](18734)'],
    products = ['[CH]=C=CO[CH][C]=C(21927)'],
    transitionState = 'TS11',
    kinetics = Arrhenius(A=(81155.9,'m^3/(mol*s)'), n=0.702818, Ea=(0,'kJ/mol'), T0=(1,'K'), comment="""Estimated using an average for rate rule [Cd_allenic;Y_rad]
Euclidian distance = 0
Multiplied by reaction path degeneracy 2.0
family: R_Recombination
Ea raised from -11.4 to 0 kJ/mol."""),
)

reaction(
    label = 'reaction12',
    reactants = ['H(8)', '[CH]=C=[C]O[CH][C]=C(25028)'],
    products = ['[CH]=C=CO[CH][C]=C(21927)'],
    transitionState = 'TS12',
    kinetics = Arrhenius(A=(2.048e+13,'cm^3/(mol*s)'), n=0.206, Ea=(0,'kJ/mol'), T0=(1,'K'), comment="""From training reaction 44 used for Cd_allenic;H_rad
Exact match found for rate rule [Cd_allenic;H_rad]
Euclidian distance = 0
family: R_Recombination
Ea raised from -0.7 to 0 kJ/mol."""),
)

reaction(
    label = 'reaction13',
    reactants = ['H(8)', '[CH]=[C][CH]OC=C=[CH](22637)'],
    products = ['[CH]=C=CO[CH][C]=C(21927)'],
    transitionState = 'TS13',
    kinetics = Arrhenius(A=(8.68156e+06,'m^3/(mol*s)'), n=0.278577, Ea=(0,'kJ/mol'), T0=(1,'K'), comment="""Estimated using an average for rate rule [Y_rad;H_rad]
Euclidian distance = 0
Multiplied by reaction path degeneracy 2.0
family: R_Recombination
Ea raised from -1.4 to 0 kJ/mol."""),
)

reaction(
    label = 'reaction14',
    reactants = ['[CH]=C=CO[CH][C]=C(21927)'],
    products = ['[CH2][C]=COC1[C]=C1(25029)'],
    transitionState = 'TS14',
    kinetics = Arrhenius(A=(4.2354e+12,'s^-1'), n=-0.1205, Ea=(181.528,'kJ/mol'), T0=(1,'K'), comment="""Estimated using average of templates [Rn;doublebond_intra_CdCdd;radadd_intra_cdsingleH] + [R3;doublebond_intra_CdCdd;radadd_intra_cdsingle] for rate rule [R3;doublebond_intra_CdCdd;radadd_intra_cdsingleH]
Euclidian distance = 1.0
family: Intra_R_Add_Endocyclic
Ea raised from 181.0 to 181.5 kJ/mol to match endothermicity of reaction."""),
)

reaction(
    label = 'reaction15',
    reactants = ['[CH]=C=CO[CH][C]=C(21927)'],
    products = ['[CH]C1=COC1[C]=C(24944)'],
    transitionState = 'TS15',
    kinetics = Arrhenius(A=(8.33681e+07,'s^-1'), n=1.1093, Ea=(126.639,'kJ/mol'), T0=(1,'K'), comment="""Estimated using an average for rate rule [R4_S_D;doublebond_intra;radadd_intra_cs]
Euclidian distance = 0
family: Intra_R_Add_Endocyclic"""),
)

reaction(
    label = 'reaction16',
    reactants = ['[CH]=C=CO[CH][C]=C(21927)'],
    products = ['[CH]C1=COC=C1[CH2](24862)'],
    transitionState = 'TS16',
    kinetics = Arrhenius(A=(6.76303e+09,'s^-1'), n=0.532663, Ea=(43.317,'kJ/mol'), T0=(1,'K'), comment="""Estimated using an average for rate rule [R5_linear;doublebond_intra;radadd_intra_cddouble]
Euclidian distance = 0
family: Intra_R_Add_Endocyclic"""),
)

reaction(
    label = 'reaction17',
    reactants = ['[CH]=C=C[O](8556)', '[CH][C]=C(18825)'],
    products = ['[CH]=C=CO[CH][C]=C(21927)'],
    transitionState = 'TS17',
    kinetics = Arrhenius(A=(43.5839,'m^3/(mol*s)'), n=1.88017, Ea=(5.1666,'kJ/mol'), T0=(1,'K'), comment="""Estimated using template [O_sec_rad;Birad] for rate rule [O_rad/OneDe;Birad]
Euclidian distance = 1.0
family: Birad_R_Recombination"""),
)

reaction(
    label = 'reaction18',
    reactants = ['H(8)', '[CH]=[C][CH]OC#C[CH2](25030)'],
    products = ['[CH]=C=CO[CH][C]=C(21927)'],
    transitionState = 'TS18',
    kinetics = Arrhenius(A=(1e+07,'m^3/(mol*s)'), n=0, Ea=(0,'kJ/mol'), T0=(1,'K'), comment="""Estimated using an average for rate rule [H_rad;Birad]
Euclidian distance = 0
family: Birad_R_Recombination"""),
)

reaction(
    label = 'reaction19',
    reactants = ['H(8)', '[C]=C=CO[CH][C]=C(25031)'],
    products = ['[CH]=C=CO[CH][C]=C(21927)'],
    transitionState = 'TS19',
    kinetics = Arrhenius(A=(1e+07,'m^3/(mol*s)'), n=0, Ea=(0,'kJ/mol'), T0=(1,'K'), comment="""Estimated using an average for rate rule [H_rad;Birad]
Euclidian distance = 0
family: Birad_R_Recombination"""),
)

reaction(
    label = 'reaction20',
    reactants = ['[CH]=C=CO[CH][C]=C(21927)'],
    products = ['[CH]=[C]C1C[C]=CO1(25032)'],
    transitionState = 'TS20',
    kinetics = Arrhenius(A=(4.13873e+09,'s^-1'), n=0.337103, Ea=(111.427,'kJ/mol'), T0=(1,'K'), comment="""Estimated using an average for rate rule [R6;doublebond_intra;radadd_intra_cs2H]
Euclidian distance = 0
family: Intra_R_Add_Exocyclic"""),
)

reaction(
    label = 'reaction21',
    reactants = ['H(8)', '[CH]=C=COC#C[CH2](25033)'],
    products = ['[CH]=C=CO[CH][C]=C(21927)'],
    transitionState = 'TS21',
    kinetics = Arrhenius(A=(4278.27,'m^3/(mol*s)'), n=1.383, Ea=(15.1097,'kJ/mol'), T0=(1,'K'), comment="""Estimated using template [Ct_Ct;HJ] for rate rule [Ct-O_Ct;HJ]
Euclidian distance = 1.0
family: R_Addition_MultipleBond"""),
)

reaction(
    label = 'reaction22',
    reactants = ['[CH]=C=C[O](8556)', 'C#C[CH2](17441)'],
    products = ['[CH]=C=CO[CH][C]=C(21927)'],
    transitionState = 'TS22',
    kinetics = Arrhenius(A=(260000,'m^3/(mol*s)'), n=0, Ea=(46.5868,'kJ/mol'), T0=(1,'K'), comment="""Estimated using average of templates [R_R;O_rad/OneDe] + [Ct_Ct;OJ_sec] for rate rule [Ct_Ct;O_rad/OneDe]
Euclidian distance = 2.0
family: R_Addition_MultipleBond"""),
)

reaction(
    label = 'reaction23',
    reactants = ['[CH]=C=CO[C]=C[CH2](25034)'],
    products = ['[CH]=C=CO[CH][C]=C(21927)'],
    transitionState = 'TS23',
    kinetics = Arrhenius(A=(1.43441e+10,'s^-1'), n=1.037, Ea=(204.137,'kJ/mol'), T0=(1,'K'), comment="""Estimated using an average for rate rule [R2H_D;Cd_rad_out_singleNd;XH_out]
Euclidian distance = 0
family: intra_H_migration"""),
)

reaction(
    label = 'reaction24',
    reactants = ['[CH]=[C][CH]OC#CC(25035)'],
    products = ['[CH]=C=CO[CH][C]=C(21927)'],
    transitionState = 'TS24',
    kinetics = Arrhenius(A=(7.74e+09,'s^-1'), n=1.08, Ea=(161.921,'kJ/mol'), T0=(1,'K'), comment="""Estimated using template [R3Hall;Cd_rad_out_singleNd;Cs_H_out_2H] for rate rule [R3HJ;Cd_rad_out_singleNd;Cs_H_out_2H]
Euclidian distance = 1.0
Multiplied by reaction path degeneracy 3.0
family: intra_H_migration"""),
)

reaction(
    label = 'reaction25',
    reactants = ['[CH]=C=[C]OC=[C]C(25036)'],
    products = ['[CH]=C=CO[CH][C]=C(21927)'],
    transitionState = 'TS25',
    kinetics = Arrhenius(A=(408000,'s^-1'), n=1.9199, Ea=(33.0402,'kJ/mol'), T0=(1,'K'), comment="""Estimated using template [R5Hall;Cd_rad_out;Cs_H_out_2H] for rate rule [R5HJ_3;Cd_rad_out_double;Cs_H_out_2H]
Euclidian distance = 1.41421356237
Multiplied by reaction path degeneracy 3.0
family: intra_H_migration"""),
)

reaction(
    label = 'reaction26',
    reactants = ['[CH]=C=CO[CH][C]=C(21927)'],
    products = ['[CH]=C1[CH]OC=[C]C1(24777)'],
    transitionState = 'TS26',
    kinetics = Arrhenius(A=(3.26172e+08,'s^-1'), n=0.58655, Ea=(80.0039,'kJ/mol'), T0=(1,'K'), comment="""Estimated using an average for rate rule [R6_linear;doublebond_intra;radadd_intra_cs2H]
Euclidian distance = 0
family: Intra_R_Add_Endocyclic"""),
)

reaction(
    label = 'reaction27',
    reactants = ['[CH]=[C]CO[C]=[C][CH2](25037)'],
    products = ['[CH]=C=CO[CH][C]=C(21927)'],
    transitionState = 'TS27',
    kinetics = Arrhenius(A=(1.4874e+09,'s^-1'), n=1.045, Ea=(63.4002,'kJ/mol'), T0=(1,'K'), comment="""Estimated using an average for rate rule [R3radExo;Y_rad;XH_Rrad]
Euclidian distance = 0
Multiplied by reaction path degeneracy 2.0
family: Intra_Disproportionation"""),
)

reaction(
    label = 'reaction28',
    reactants = ['[CH]=C[CH]O[C]=[C][CH2](25038)'],
    products = ['[CH]=C=CO[CH][C]=C(21927)'],
    transitionState = 'TS28',
    kinetics = Arrhenius(A=(5.14222e+08,'s^-1'), n=0.311, Ea=(24.9733,'kJ/mol'), T0=(1,'K'), comment="""Estimated using template [R4;Y_rad;XH_Rrad] for rate rule [R4radEndo;Y_rad;XH_Rrad]
Euclidian distance = 1.0
family: Intra_Disproportionation"""),
)

reaction(
    label = 'reaction49',
    reactants = ['[CH]=C=CO[CH][C]=C(21927)'],
    products = ['[CH]C(C=O)=C[C]=C(22615)'],
    transitionState = 'TS29',
    kinetics = Arrhenius(A=(7040,'s^-1'), n=2.66, Ea=(313.8,'kJ/mol'), T0=(1,'K'), Tmin=(300,'K'), Tmax=(1500,'K'), comment="""Estimated using an average for rate rule [R_ROR;R1_doublebond;R2_doublebond_H;R_O_C]
Euclidian distance = 0
family: ketoenol"""),
)

reaction(
    label = 'reaction30',
    reactants = ['C#CC=O(21959)', '[CH][C]=C(18825)'],
    products = ['[CH]=C=CO[CH][C]=C(21927)'],
    transitionState = 'TS30',
    kinetics = Arrhenius(A=(1.07e+06,'cm^3/(mol*s)'), n=2.43, Ea=(22.5936,'kJ/mol'), T0=(1,'K'), Tmin=(300,'K'), Tmax=(1500,'K'), comment="""Estimated using template [Od_CO-CtH;YJ] for rate rule [Od_CO-CtH;CJ]
Euclidian distance = 1.0
family: R_Addition_MultipleBond"""),
)

reaction(
    label = 'reaction31',
    reactants = ['[CH]=[C]COC#C[CH2](24817)'],
    products = ['[CH]=C=CO[CH][C]=C(21927)'],
    transitionState = 'TS31',
    kinetics = Arrhenius(A=(3.32e+07,'s^-1'), n=1.69, Ea=(159.41,'kJ/mol'), T0=(1,'K'), comment="""Estimated using template [R3H_SS_O;Y_rad_out;Cs_H_out_H/OneDe] for rate rule [R3H_SS_O;Cd_rad_out;Cs_H_out_H/Ct]
Euclidian distance = 1.41421356237
Multiplied by reaction path degeneracy 2.0
family: intra_H_migration"""),
)

reaction(
    label = 'reaction32',
    reactants = ['[C]#CCOC=[C][CH2](25039)'],
    products = ['[CH]=C=CO[CH][C]=C(21927)'],
    transitionState = 'TS32',
    kinetics = Arrhenius(A=(1.47715e+10,'s^-1'), n=0.8, Ea=(147.277,'kJ/mol'), T0=(1,'K'), comment="""Estimated using template [R3H_MS;Y_rad_out;Cs_H_out_1H] for rate rule [R3H_TS;Ct_rad_out;Cs_H_out_H/NonDeO]
Euclidian distance = 1.73205080757
Multiplied by reaction path degeneracy 2.0
family: intra_H_migration"""),
)

reaction(
    label = 'reaction33',
    reactants = ['[C]=C=CO[CH]C=C(25040)'],
    products = ['[CH]=C=CO[CH][C]=C(21927)'],
    transitionState = 'TS33',
    kinetics = Arrhenius(A=(19101.7,'s^-1'), n=1.79759, Ea=(74.4436,'kJ/mol'), T0=(1,'K'), comment="""Estimated using template [R6Hall;Y_rad_out;XH_out] for rate rule [R6HJ_2;Ct_rad_out;XH_out]
Euclidian distance = 1.41421356237
family: intra_H_migration"""),
)

reaction(
    label = 'reaction34',
    reactants = ['[C]#C[CH]OC=[C]C(25041)'],
    products = ['[CH]=C=CO[CH][C]=C(21927)'],
    transitionState = 'TS34',
    kinetics = Arrhenius(A=(305745,'s^-1'), n=1.095, Ea=(83.5753,'kJ/mol'), T0=(1,'K'), comment="""Estimated using template [R7Hall;Y_rad_out;Cs_H_out_2H] for rate rule [R7Hall;Ct_rad_out;Cs_H_out_2H]
Euclidian distance = 1.0
Multiplied by reaction path degeneracy 3.0
family: intra_H_migration"""),
)

reaction(
    label = 'reaction35',
    reactants = ['[CH]=C=CO[CH][C]=C(21927)'],
    products = ['[CH2]C1[CH]OC=[C]C=1(24893)'],
    transitionState = 'TS35',
    kinetics = Arrhenius(A=(9.23539e+09,'s^-1'), n=0.445806, Ea=(31.0324,'kJ/mol'), T0=(1,'K'), comment="""Estimated using an average for rate rule [R6_linear;triplebond_intra_H;radadd_intra]
Euclidian distance = 0
family: Intra_R_Add_Endocyclic"""),
)

reaction(
    label = 'reaction36',
    reactants = ['[CH]=C=CO[CH][C]=C(21927)'],
    products = ['[C]1[CH]OC=[C]CC=1(24811)'],
    transitionState = 'TS36',
    kinetics = Arrhenius(A=(1.8036e+08,'s^-1'), n=0.568448, Ea=(107.224,'kJ/mol'), T0=(1,'K'), comment="""Estimated using average of templates [R7_linear;multiplebond_intra;radadd_intra_cs2H] + [R7_linear;triplebond_intra_H;radadd_intra] for rate rule [R7_linear;triplebond_intra_H;radadd_intra_cs2H]
Euclidian distance = 2.0
family: Intra_R_Add_Endocyclic"""),
)

reaction(
    label = 'reaction37',
    reactants = ['[CH]=C=CO[CH][C]=C(21927)'],
    products = ['C#CC1OC=C1[CH2](24933)'],
    transitionState = 'TS37',
    kinetics = Arrhenius(A=(2e+12,'s^-1'), n=0, Ea=(7.5312,'kJ/mol'), T0=(1,'K'), comment="""Estimated using template [Rn;C_rad_out_H/OneDe;Ypri_rad_out] for rate rule [R4_SSD;C_rad_out_H/OneDe;Ypri_rad_out]
Euclidian distance = 2.0
family: Birad_recombination"""),
)

reaction(
    label = 'reaction38',
    reactants = ['[C]#C(5143)', '[CH]OC=[C][CH2](19096)'],
    products = ['[CH]=C=CO[CH][C]=C(21927)'],
    transitionState = 'TS38',
    kinetics = Arrhenius(A=(1.14854e+06,'m^3/(mol*s)'), n=0.575199, Ea=(34.3157,'kJ/mol'), T0=(1,'K'), comment="""Estimated using template [Y_rad;Birad] for rate rule [Ct_rad;Birad]
Euclidian distance = 1.0
family: Birad_R_Recombination"""),
)

reaction(
    label = 'reaction39',
    reactants = ['[C]=C=COC[C]=C(24330)'],
    products = ['[CH]=C=CO[CH][C]=C(21927)'],
    transitionState = 'TS39',
    kinetics = Arrhenius(A=(2.20409e+06,'s^-1'), n=1.54368, Ea=(52.1315,'kJ/mol'), T0=(1,'K'), comment="""Estimated using template [R5Hall;Y_rad_out;XH_out] for rate rule [R5HJ_2;Ct_rad_out;XH_out]
Euclidian distance = 1.41421356237
Multiplied by reaction path degeneracy 2.0
family: intra_H_migration"""),
)

reaction(
    label = 'reaction40',
    reactants = ['[CH]=C=CO[CH][C]=C(21927)'],
    products = ['C=[C]C1[CH][C]=CO1(24991)'],
    transitionState = 'TS40',
    kinetics = Arrhenius(A=(1.21395e+11,'s^-1'), n=0.177667, Ea=(55.8076,'kJ/mol'), T0=(1,'K'), comment="""Estimated using template [R5_linear;triplebond_intra_H;radadd_intra] for rate rule [R5_linear;triplebond_intra_H;radadd_intra_cs]
Euclidian distance = 1.0
family: Intra_R_Add_Endocyclic"""),
)

reaction(
    label = 'reaction41',
    reactants = ['[CH]=C=CO[CH][C]=C(21927)'],
    products = ['C#CC1OC1[C]=C(25042)'],
    transitionState = 'TS41',
    kinetics = Arrhenius(A=(1.61967e+11,'s^-1'), n=0.0247333, Ea=(7.86034,'kJ/mol'), T0=(1,'K'), comment="""Estimated using template [Rn;Y_rad_out;Cpri_rad_out_single] for rate rule [R3_SS;Y_rad_out;Cpri_rad_out_H/OneDe]
Euclidian distance = 2.82842712475
family: Birad_recombination"""),
)

reaction(
    label = 'reaction42',
    reactants = ['[CH]=C[CH]OC#C[CH2](25043)'],
    products = ['[CH]=C=CO[CH][C]=C(21927)'],
    transitionState = 'TS42',
    kinetics = Arrhenius(A=(1.13341e+09,'s^-1'), n=1.04717, Ea=(150.072,'kJ/mol'), T0=(1,'K'), comment="""Estimated using template [R4Hall;Cd_rad_out_double;XH_out] for rate rule [R4HJ_2;Cd_rad_out_double;XH_out]
Euclidian distance = 1.0
family: intra_H_migration"""),
)

reaction(
    label = 'reaction43',
    reactants = ['[CH][C]=C[O](21209)', 'C#C[CH2](17441)'],
    products = ['[CH]=C=CO[CH][C]=C(21927)'],
    transitionState = 'TS43',
    kinetics = Arrhenius(A=(40577.9,'m^3/(mol*s)'), n=0.702818, Ea=(0,'kJ/mol'), T0=(1,'K'), comment="""Estimated using an average for rate rule [Y_rad;Cd_allenic]
Euclidian distance = 0
family: R_Recombination
Ea raised from -11.4 to 0 kJ/mol."""),
)

reaction(
    label = 'reaction44',
    reactants = ['[CH]=[C][CH]O[C]=C[CH2](21941)'],
    products = ['[CH]=C=CO[CH][C]=C(21927)'],
    transitionState = 'TS44',
    kinetics = Arrhenius(A=(2.24409e+10,'s^-1'), n=0.34095, Ea=(22.3009,'kJ/mol'), T0=(1,'K'), comment="""Estimated using average of templates [Rn;Y_rad_NDe;XH_Rrad] + [R2radExo;Y_rad;XH_Rrad] for rate rule [R2radExo;Y_rad_NDe;XH_Rrad]
Euclidian distance = 1.0
family: Intra_Disproportionation"""),
)

reaction(
    label = 'reaction45',
    reactants = ['[CH]C=CO[CH][C]=[CH](21936)'],
    products = ['[CH]=C=CO[CH][C]=C(21927)'],
    transitionState = 'TS45',
    kinetics = Arrhenius(A=(1.949e+11,'s^-1'), n=0.486, Ea=(22.8614,'kJ/mol'), T0=(1,'K'), comment="""Estimated using an average for rate rule [R2radExo;Y_rad;XH_Rrad_NDe]
Euclidian distance = 0
family: Intra_Disproportionation"""),
)

reaction(
    label = 'reaction46',
    reactants = ['C=C=C[O](12570)', '[CH][C]=[CH](21256)'],
    products = ['[CH]=C=CO[CH][C]=C(21927)'],
    transitionState = 'TS46',
    kinetics = Arrhenius(A=(87.1677,'m^3/(mol*s)'), n=1.88017, Ea=(5.1666,'kJ/mol'), T0=(1,'K'), comment="""Estimated using template [O_sec_rad;Birad] for rate rule [O_rad/OneDe;Birad]
Euclidian distance = 1.0
Multiplied by reaction path degeneracy 2.0
family: Birad_R_Recombination"""),
)

reaction(
    label = 'reaction47',
    reactants = ['H(8)', '[CH]=C=[C]OC=C=C(25044)'],
    products = ['[CH]=C=CO[CH][C]=C(21927)'],
    transitionState = 'TS47',
    kinetics = Arrhenius(A=(4278.27,'m^3/(mol*s)'), n=1.383, Ea=(15.1097,'kJ/mol'), T0=(1,'K'), comment="""Estimated using template [Ct_Ct;HJ] for rate rule [Ct-O_Ct;HJ]
Euclidian distance = 1.0
family: R_Addition_MultipleBond"""),
)

reaction(
    label = 'reaction48',
    reactants = ['[CH]C=[C]OC=C=C(25045)'],
    products = ['[CH]=C=CO[CH][C]=C(21927)'],
    transitionState = 'TS48',
    kinetics = Arrhenius(A=(1.43441e+10,'s^-1'), n=1.037, Ea=(204.137,'kJ/mol'), T0=(1,'K'), comment="""Estimated using an average for rate rule [R2H_D;Cd_rad_out_singleNd;XH_out]
Euclidian distance = 0
family: intra_H_migration"""),
)

reaction(
    label = 'reaction49',
    reactants = ['[CH][C]=[C]OC[C]=C(25046)'],
    products = ['[CH]=C=CO[CH][C]=C(21927)'],
    transitionState = 'TS49',
    kinetics = Arrhenius(A=(4.00798e+09,'s^-1'), n=0.37, Ea=(78.2471,'kJ/mol'), T0=(1,'K'), comment="""Estimated using average of templates [R3;Y_rad;XH_Rrad_De] + [R3radExo;Y_rad;XH_Rrad] for rate rule [R3radExo;Y_rad;XH_Rrad_De]
Euclidian distance = 1.0
Multiplied by reaction path degeneracy 2.0
family: Intra_Disproportionation"""),
)

reaction(
    label = 'reaction50',
    reactants = ['[CH][C]=[C]O[CH]C=C(21935)'],
    products = ['[CH]=C=CO[CH][C]=C(21927)'],
    transitionState = 'TS50',
    kinetics = Arrhenius(A=(5.14222e+08,'s^-1'), n=0.311, Ea=(24.9733,'kJ/mol'), T0=(1,'K'), comment="""Estimated using template [R4;Y_rad;XH_Rrad] for rate rule [R4radEndo;Y_rad;XH_Rrad]
Euclidian distance = 1.0
family: Intra_Disproportionation"""),
)

reaction(
    label = 'reaction51',
    reactants = ['[CH]=C=[C]O[CH][C]C(25047)'],
    products = ['[CH]=C=CO[CH][C]=C(21927)'],
    transitionState = 'TS51',
    kinetics = Arrhenius(A=(6.37831e+09,'s^-1'), n=0.137, Ea=(24.9733,'kJ/mol'), T0=(1,'K'), comment="""Estimated using template [R5;Y_rad;XH_Rrad] for rate rule [R5radEndo;Y_rad;XH_Rrad]
Euclidian distance = 1.0
Multiplied by reaction path degeneracy 3.0
family: Intra_Disproportionation"""),
)

reaction(
    label = 'reaction52',
    reactants = ['[CH]=C=CO[CH][C]=C(21927)'],
    products = ['[CH]=C=CC([CH2])=C[O](21884)'],
    transitionState = 'TS52',
    kinetics = Arrhenius(A=(7040,'s^-1'), n=2.66, Ea=(313.8,'kJ/mol'), T0=(1,'K'), Tmin=(300,'K'), Tmax=(1500,'K'), comment="""Estimated using an average for rate rule [R_ROR;R1_doublebond;R2_doublebond_H;R_O_C]
Euclidian distance = 0
family: ketoenol"""),
)

network(
    label = '4916',
    isomers = [
        '[CH]=C=CO[CH][C]=C(21927)',
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
    label = '4916',
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

