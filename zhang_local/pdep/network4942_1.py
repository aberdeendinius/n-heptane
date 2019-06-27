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
    label = '[CH]C1=COC1[C]=[CH](25320)',
    structure = SMILES('[CH]C1=COC1[C]=[CH]'),
    E0 = (879.316,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
    molecularWeight = (92.0954,'amu'),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[0.98528,0.0466654,1.65101e-05,-6.87926e-08,3.34374e-11,105884,23.8959], Tmin=(100,'K'), Tmax=(937.498,'K')), NASAPolynomial(coeffs=[21.7085,0.00791719,-9.66047e-07,1.49492e-10,-1.77675e-14,99815.4,-86.3856], Tmin=(937.498,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(879.316,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(249.434,'J/(mol*K)'), comment="""Thermo group additivity estimation: group(O2s-Cs(Cds-Cd)) + group(Cs-(Cds-Cds)(Cds-Cds)OsH) + group(Cs-(Cds-Cds)HHH) + group(Cds-CdsCsCs) + group(Cds-CdsCsH) + group(Cds-CdsOsH) + group(Cds-CdsHH) + ring(Cyclobutene) + radical(Cds_P) + radical(Cds_S) + radical(AllylJ2_triplet)"""),
)

species(
    label = '[CH]C1([CH][O])C=C=C1(26121)',
    structure = SMILES('[CH]C1([CH][O])C=C=C1'),
    E0 = (1134.69,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
    molecularWeight = (92.0954,'amu'),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[1.34377,0.0653323,-9.89629e-05,9.15168e-08,-3.46487e-11,136560,23.5624], Tmin=(100,'K'), Tmax=(764.228,'K')), NASAPolynomial(coeffs=[5.4074,0.0340277,-1.78221e-05,3.55185e-09,-2.52053e-13,136232,6.9677], Tmin=(764.228,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(1134.69,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(249.434,'J/(mol*K)'), comment="""Thermo group additivity estimation: group(O2s-CsH) + group(Cs-(Cds-Cds)(Cds-Cds)CsCs) + group(Cs-CsOsHH) + group(Cs-CsHHH) + group(Cds-CdsCsH) + group(Cds-CdsCsH) + group(Cdd-CdsCds) + ring(cyclobutadiene_13) + radical(CCJ2_triplet) + radical(CCOJ) + radical(CCsJOH)"""),
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
    label = '[CH]C([C]=O)=CC#C(26122)',
    structure = SMILES('[CH]C([C]=O)=CC#C'),
    E0 = (638.164,'kJ/mol'),
    modes = [
        HarmonicOscillator(frequencies=([1855,455,950,750,770,3400,2100,350,440,435,1725,2175,525,3010,987.5,1337.5,450,1655,180,180,180],'cm^-1')),
        HinderedRotor(inertia=(2.04402,'amu*angstrom^2'), symmetry=1, barrier=(46.9959,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(2.04454,'amu*angstrom^2'), symmetry=1, barrier=(47.0079,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(2.04466,'amu*angstrom^2'), symmetry=1, barrier=(47.0109,'kJ/mol'), semiclassical=False),
    ],
    spinMultiplicity = 4,
    opticalIsomers = 1,
    molecularWeight = (91.0874,'amu'),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[1.57951,0.0572997,-5.29201e-05,2.57249e-08,-5.28758e-12,76837.1,21.6861], Tmin=(100,'K'), Tmax=(1124.63,'K')), NASAPolynomial(coeffs=[9.70993,0.0283824,-1.43513e-05,2.86201e-09,-2.0533e-13,75008.4,-18.4935], Tmin=(1124.63,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(638.164,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(220.334,'J/(mol*K)'), comment="""Thermo group additivity estimation: group(Cs-(Cds-Cds)HHH) + group(Cd-CdCs(CO)) + group(Cds-CdsCtH) + group(Cds-O2d(Cds-Cds)H) + group(Ct-Ct(Cds-Cds)) + group(Ct-CtH) + radical(C=C(C)CJ=O) + radical(AllylJ2_triplet)"""),
)

species(
    label = '[CH]C(=[C]O)C=C=[CH](26123)',
    structure = SMILES('[CH]C(=[C]O)C=C=[CH]'),
    E0 = (753.022,'kJ/mol'),
    modes = [
        HarmonicOscillator(frequencies=([3120,650,792.5,1650,3615,1277.5,1000,540,610,2055,350,440,435,1725,1685,370,3010,987.5,1337.5,450,1655,180,180,180],'cm^-1')),
        HinderedRotor(inertia=(2.02892,'amu*angstrom^2'), symmetry=1, barrier=(46.6489,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(2.02491,'amu*angstrom^2'), symmetry=1, barrier=(46.5566,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(2.0268,'amu*angstrom^2'), symmetry=1, barrier=(46.6,'kJ/mol'), semiclassical=False),
    ],
    spinMultiplicity = 5,
    opticalIsomers = 1,
    molecularWeight = (92.0954,'amu'),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[-0.170183,0.077595,-8.65942e-05,4.69947e-08,-9.52668e-12,90729.9,29.1404], Tmin=(100,'K'), Tmax=(1389.74,'K')), NASAPolynomial(coeffs=[19.5748,0.00934815,-6.10878e-07,-1.62867e-10,1.97821e-14,86344.3,-68.6498], Tmin=(1389.74,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(753.022,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(245.277,'J/(mol*K)'), comment="""Thermo group additivity estimation: group(O2s-(Cds-Cd)H) + group(Cs-(Cds-Cds)HHH) + group(Cds-Cds(Cds-Cds)Cs) + group(Cds-Cds(Cds-Cds)H) + group(Cds-CdsOsH) + group(Cds-CdsHH) + group(Cdd-CdsCds) + radical(AllylJ2_triplet) + radical(C=CJO) + radical(C=C=CJ)"""),
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
    label = '[CH]C(=C=C=[CH])[CH]O(26124)',
    structure = SMILES('[CH]C(=C=C=[CH])[CH]O'),
    E0 = (737.547,'kJ/mol'),
    modes = [
        HarmonicOscillator(frequencies=([3120,650,792.5,1650,3025,407.5,1350,352.5,540,563.333,586.667,610,1970,2140,3615,1277.5,1000,350,440,435,1725,180,180,180],'cm^-1')),
        HinderedRotor(inertia=(2.16463,'amu*angstrom^2'), symmetry=1, barrier=(49.769,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(2.16456,'amu*angstrom^2'), symmetry=1, barrier=(49.7675,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(2.16586,'amu*angstrom^2'), symmetry=1, barrier=(49.7974,'kJ/mol'), semiclassical=False),
    ],
    spinMultiplicity = 5,
    opticalIsomers = 1,
    molecularWeight = (92.0954,'amu'),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[0.809817,0.0644828,-5.90132e-05,2.57244e-08,-3.63864e-12,88826.2,27.1788], Tmin=(100,'K'), Tmax=(987.564,'K')), NASAPolynomial(coeffs=[14.2468,0.0202388,-7.27478e-06,1.23645e-09,-8.20786e-14,85675.8,-39.9921], Tmin=(987.564,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(737.547,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(245.277,'J/(mol*K)'), comment="""Thermo group additivity estimation: group(O2s-CsH) + group(Cs-(Cds-Cds)OsHH) + group(Cs-(Cds-Cds)HHH) + group(Cds-CdsCsCs) + group(Cds-CdsHH) + group(Cdd-CdsCds) + group(Cdd-CdsCds) + radical(AllylJ2_triplet) + radical(C=CCJO) + radical(C=C=CJ)"""),
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
    label = '[CH]C([C]=C=[CH])=C[O](26125)',
    structure = SMILES('[CH]C([C]=C=[CH])=C[O]'),
    E0 = (892.582,'kJ/mol'),
    modes = [
        HarmonicOscillator(frequencies=([3120,650,792.5,1650,1685,370,540,610,2055,350,440,435,1725,3010,987.5,1337.5,450,1655,180,180,180,180],'cm^-1')),
        HinderedRotor(inertia=(2.03708,'amu*angstrom^2'), symmetry=1, barrier=(46.8364,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(2.03554,'amu*angstrom^2'), symmetry=1, barrier=(46.8012,'kJ/mol'), semiclassical=False),
    ],
    spinMultiplicity = 6,
    opticalIsomers = 1,
    molecularWeight = (91.0874,'amu'),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[0.0022421,0.073372,-8.05368e-05,4.2658e-08,-8.44187e-12,107509,27.2801], Tmin=(100,'K'), Tmax=(1420.09,'K')), NASAPolynomial(coeffs=[19.3295,0.00832222,-6.19464e-07,-1.20726e-10,1.53469e-14,103090,-68.9746], Tmin=(1420.09,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(892.582,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(224.491,'J/(mol*K)'), comment="""Thermo group additivity estimation: group(O2s-(Cds-Cd)H) + group(Cs-(Cds-Cds)HHH) + group(Cds-Cds(Cds-Cds)Cs) + group(Cds-Cds(Cds-Cds)H) + group(Cds-CdsOsH) + group(Cds-CdsHH) + group(Cdd-CdsCds) + radical(C=C=CJ) + radical(Cds_S) + radical(C=COJ) + radical(AllylJ2_triplet)"""),
)

species(
    label = '[CH][C]([C]=O)C=C=[CH](26126)',
    structure = SMILES('[CH][C]([C]=O)C=C=[CH]'),
    E0 = (891.663,'kJ/mol'),
    modes = [
        HarmonicOscillator(frequencies=([1855,455,950,360,370,350,3010,987.5,1337.5,450,1655,3120,650,792.5,1650,540,610,2055,274.073,783.405,783.464],'cm^-1')),
        HinderedRotor(inertia=(0.0181267,'amu*angstrom^2'), symmetry=1, barrier=(82.2003,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.188647,'amu*angstrom^2'), symmetry=1, barrier=(82.1971,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(3.57518,'amu*angstrom^2'), symmetry=1, barrier=(82.2005,'kJ/mol'), semiclassical=False),
    ],
    spinMultiplicity = 6,
    opticalIsomers = 1,
    molecularWeight = (91.0874,'amu'),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[1.34372,0.0523784,-5.12675e-05,2.57282e-08,-5.07191e-12,107343,26.189], Tmin=(100,'K'), Tmax=(1240.97,'K')), NASAPolynomial(coeffs=[13.2452,0.0140158,-4.8965e-06,8.16588e-10,-5.32299e-14,104389,-33.7979], Tmin=(1240.97,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(891.663,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(220.334,'J/(mol*K)'), comment="""Thermo group additivity estimation: group(Cs-(Cds-O2d)CsCsH) + group(Cs-CsHHH) + group(Cds-CdsCsH) + group(Cds-OdCsH) + group(Cds-CdsHH) + group(Cdd-CdsCds) + radical(C=C=CJ) + radical(CC(C)CJ=O) + radical(CCJ2_triplet) + radical(C=CCJ(C)C=O)"""),
)

species(
    label = '[CH]C1([CH]O1)C=C=[CH](26127)',
    structure = SMILES('[CH]C1([CH]O1)C=C=[CH]'),
    E0 = (878.408,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
    molecularWeight = (92.0954,'amu'),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[-0.312022,0.0736567,-8.60784e-05,4.81123e-08,-9.62738e-12,105822,29.0654], Tmin=(100,'K'), Tmax=(1519.15,'K')), NASAPolynomial(coeffs=[17.2851,0.00481307,4.12317e-06,-1.22541e-09,9.70085e-14,103073,-54.6397], Tmin=(1519.15,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(878.408,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(249.434,'J/(mol*K)'), comment="""Thermo group additivity estimation: group(O2s-CsCs) + group(Cs-(Cds-Cds)CsCsOs) + group(Cs-CsOsHH) + group(Cs-CsHHH) + group(Cds-CdsCsH) + group(Cds-CdsHH) + group(Cdd-CdsCds) + ring(Ethylene_oxide) + radical(CCsJO) + radical(CCJ2_triplet) + radical(C=C=CJ)"""),
)

species(
    label = '[CH]C(=C[O])C1[C]=C1(26128)',
    structure = SMILES('[CH]C(=C[O])C1[C]=C1'),
    E0 = (843.375,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
    molecularWeight = (92.0954,'amu'),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[0.87315,0.0615243,-4.91595e-05,1.44733e-08,4.94969e-13,101554,23.7135], Tmin=(100,'K'), Tmax=(979.103,'K')), NASAPolynomial(coeffs=[14.9374,0.0189707,-6.80018e-06,1.17825e-09,-8.00992e-14,98085.2,-47.4898], Tmin=(979.103,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(843.375,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(249.434,'J/(mol*K)'), comment="""Thermo group additivity estimation: group(O2s-(Cds-Cd)H) + group(Cs-(Cds-Cds)(Cds-Cds)(Cds-Cds)H) + group(Cs-(Cds-Cds)HHH) + group(Cds-CdsCsCs) + group(Cds-CdsCsH) + group(Cds-CdsCsH) + group(Cds-CdsOsH) + ring(Cyclopropene) + radical(C=COJ) + radical(cyclopropenyl-vinyl) + radical(AllylJ2_triplet)"""),
)

species(
    label = '[CH]C1=COC([CH])=C1(26106)',
    structure = SMILES('[CH]C1=COC([CH])=C1'),
    E0 = (614.788,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
    molecularWeight = (92.0954,'amu'),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[1.79043,0.0365206,2.06486e-05,-4.60019e-08,1.80151e-11,74031.9,22.749], Tmin=(100,'K'), Tmax=(1031.96,'K')), NASAPolynomial(coeffs=[10.3653,0.0315879,-1.33237e-05,2.52349e-09,-1.79409e-13,70755,-26.1919], Tmin=(1031.96,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(614.788,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(249.434,'J/(mol*K)'), comment="""Thermo group additivity estimation: group(O2s-(Cds-Cd)(Cds-Cd)) + group(Cs-(Cds-Cds)HHH) + group(Cs-(Cds-Cds)HHH) + group(Cds-Cds(Cds-Cds)Cs) + group(Cds-CdsCsOs) + group(Cds-Cds(Cds-Cds)H) + group(Cds-CdsOsH) + ring(Furan) + radical(AllylJ2_triplet) + radical(AllylJ2_triplet)"""),
)

species(
    label = '[CH]=C1C=[C][CH]C1[O](26129)',
    structure = SMILES('[CH]=C1C=[C][CH]C1[O]'),
    E0 = (729.257,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
    molecularWeight = (92.0954,'amu'),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[1.95086,0.0215923,6.97311e-05,-1.15627e-07,4.80286e-11,87804.6,22.6476], Tmin=(100,'K'), Tmax=(947.142,'K')), NASAPolynomial(coeffs=[19.5579,0.00665412,-7.15875e-07,1.96187e-10,-2.67581e-14,81804.1,-75.41], Tmin=(947.142,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(729.257,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(257.749,'J/(mol*K)'), comment="""Thermo group additivity estimation: group(O2s-CsH) + group(Cs-(Cds-Cds)CsOsH) + group(Cs-(Cds-Cds)CsHH) + group(Cds-Cds(Cds-Cds)Cs) + group(Cds-CdsCsH) + group(Cds-Cds(Cds-Cds)H) + group(Cds-CdsHH) + ring(3-Methylenecyclopentene) + radical(C=CCJC(O)C=C) + radical(CC(C)OJ) + radical(Cds_P) + radical(cyclopentene-vinyl)"""),
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
    label = '[CH]C([C]=C=[CH])[CH][O](26130)',
    structure = SMILES('[CH]C([C]=C=[CH])[CH][O]'),
    E0 = (1179.01,'kJ/mol'),
    modes = [
        HarmonicOscillator(frequencies=([3120,650,792.5,1650,3025,407.5,1350,352.5,540,610,2055,1380,1390,370,380,2900,435,1685,370,180,180,833.806,834.127,834.291],'cm^-1')),
        HinderedRotor(inertia=(0.161799,'amu*angstrom^2'), symmetry=1, barrier=(3.72008,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.161819,'amu*angstrom^2'), symmetry=1, barrier=(3.72053,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.161801,'amu*angstrom^2'), symmetry=1, barrier=(3.72013,'kJ/mol'), semiclassical=False),
    ],
    spinMultiplicity = 7,
    opticalIsomers = 1,
    molecularWeight = (92.0954,'amu'),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[0.717182,0.0817192,-0.000141662,1.26551e-07,-4.27293e-11,141912,29.9378], Tmin=(100,'K'), Tmax=(892.227,'K')), NASAPolynomial(coeffs=[8.0673,0.0267592,-1.22643e-05,2.22011e-09,-1.45647e-13,141476,0.224947], Tmin=(892.227,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(1179.01,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(245.277,'J/(mol*K)'), comment="""Thermo group additivity estimation: group(O2s-CsH) + group(Cs-(Cds-Cdd-Cd)CsCsH) + group(Cs-CsOsHH) + group(Cs-CsHHH) + group(Cds-CdsCsH) + group(Cds-CdsHH) + group(Cdd-CdsCds) + radical(CCJ2_triplet) + radical(CCsJOH) + radical(CCOJ) + radical(Cds_S) + radical(C=C=CJ)"""),
)

species(
    label = '[CH]C(=[C][O])C[C]=[CH](26131)',
    structure = SMILES('[CH]C(=[C][O])C[C]=[CH]'),
    E0 = (1077.73,'kJ/mol'),
    modes = [
        HarmonicOscillator(frequencies=([3120,650,792.5,1650,2750,2850,1437.5,1250,1305,750,350,350,440,435,1725,1670,1700,300,440,375.14,375.187,375.255,375.275,375.324],'cm^-1')),
        HinderedRotor(inertia=(0.512367,'amu*angstrom^2'), symmetry=1, barrier=(51.207,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.512722,'amu*angstrom^2'), symmetry=1, barrier=(51.211,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.513033,'amu*angstrom^2'), symmetry=1, barrier=(51.2088,'kJ/mol'), semiclassical=False),
    ],
    spinMultiplicity = 7,
    opticalIsomers = 1,
    molecularWeight = (92.0954,'amu'),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[1.10735,0.0638036,-6.54496e-05,3.64963e-08,-8.30762e-12,129725,30.3511], Tmin=(100,'K'), Tmax=(1055.75,'K')), NASAPolynomial(coeffs=[11.2571,0.0253474,-1.08101e-05,1.9926e-09,-1.37005e-13,127582,-19.1661], Tmin=(1055.75,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(1077.73,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(245.277,'J/(mol*K)'), comment="""Thermo group additivity estimation: group(O2s-(Cds-Cd)H) + group(Cs-(Cds-Cds)(Cds-Cds)HH) + group(Cs-(Cds-Cds)HHH) + group(Cds-CdsCsCs) + group(Cds-CdsCsH) + group(Cds-CdsOsH) + group(Cds-CdsHH) + radical(Cds_S) + radical(AllylJ2_triplet) + radical(C=CJO) + radical(Cds_P) + radical(C=COJ)"""),
)

species(
    label = '[CH]C=C[C]([CH])[C]=O(26132)',
    structure = SMILES('[CH]C=C[C]([CH])[C]=O'),
    E0 = (931.266,'kJ/mol'),
    modes = [
        HarmonicOscillator(frequencies=([360,370,350,2995,3025,975,1000,1300,1375,400,500,1630,1680,1855,455,950,200,800,960,1120,1280,1440,1600],'cm^-1')),
        HinderedRotor(inertia=(0.156089,'amu*angstrom^2'), symmetry=1, barrier=(3.5888,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.156089,'amu*angstrom^2'), symmetry=1, barrier=(3.5888,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.156089,'amu*angstrom^2'), symmetry=1, barrier=(3.5888,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.156089,'amu*angstrom^2'), symmetry=1, barrier=(3.5888,'kJ/mol'), semiclassical=False),
    ],
    spinMultiplicity = 7,
    opticalIsomers = 1,
    molecularWeight = (92.0954,'amu'),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[1.33833,0.0510147,-3.11762e-05,6.62825e-09,3.68675e-13,112108,27.643], Tmin=(100,'K'), Tmax=(1176.41,'K')), NASAPolynomial(coeffs=[11.9297,0.0248151,-1.02819e-05,1.87806e-09,-1.28668e-13,108937,-28.0609], Tmin=(1176.41,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(931.266,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(241.12,'J/(mol*K)'), comment="""Thermo group additivity estimation: group(Cs-(Cds-O2d)CsCsH) + group(Cs-CsHHH) + group(Cs-(Cds-Cds)HHH) + group(Cds-CdsCsH) + group(Cds-CdsCsH) + group(Cds-OdCsH) + radical(AllylJ2_triplet) + radical(CC(C)CJ=O) + radical(CCJ2_triplet) + radical(C=CCJ(C)C=O)"""),
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
    label = '[CH]C(=[CH])C=C=[CH](21251)',
    structure = SMILES('[CH]C([CH])=CC#C'),
    E0 = (936.996,'kJ/mol'),
    modes = [
        HarmonicOscillator(frequencies=([750,770,3400,2100,350,440,435,1725,2175,525,3010,987.5,1337.5,450,1655,180,180,180,180,180,180],'cm^-1')),
        HinderedRotor(inertia=(2.14427,'amu*angstrom^2'), symmetry=1, barrier=(49.3009,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(2.14529,'amu*angstrom^2'), symmetry=1, barrier=(49.3244,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(2.14837,'amu*angstrom^2'), symmetry=1, barrier=(49.3953,'kJ/mol'), semiclassical=False),
    ],
    spinMultiplicity = 5,
    opticalIsomers = 1,
    molecularWeight = (76.096,'amu'),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[1.26545,0.0527517,-3.08891e-05,4.8324e-09,1.317e-12,112799,21.4875], Tmin=(100,'K'), Tmax=(1100.67,'K')), NASAPolynomial(coeffs=[11.4962,0.0267993,-1.08221e-05,1.94592e-09,-1.32441e-13,109867,-31.9409], Tmin=(1100.67,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(936.996,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(220.334,'J/(mol*K)'), comment="""Thermo group additivity estimation: group(Cs-(Cds-Cds)HHH) + group(Cs-(Cds-Cds)HHH) + group(Cds-CdsCsCs) + group(Cds-CdsCtH) + group(Ct-Ct(Cds-Cds)) + group(Ct-CtH) + radical(AllylJ2_triplet) + radical(AllylJ2_triplet)"""),
)

species(
    label = '[C]#C[CH]C([CH])=C[O](26133)',
    structure = SMILES('[C]#C[CH]C([CH])=C[O]'),
    E0 = (1054.11,'kJ/mol'),
    modes = [
        HarmonicOscillator(frequencies=([3025,407.5,1350,352.5,350,440,435,1725,2175,525,3010,987.5,1337.5,450,1655,428.024,428.024,428.028,428.03,428.03,428.035],'cm^-1')),
        HinderedRotor(inertia=(0.396364,'amu*angstrom^2'), symmetry=1, barrier=(51.5293,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.396363,'amu*angstrom^2'), symmetry=1, barrier=(51.5294,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.396354,'amu*angstrom^2'), symmetry=1, barrier=(51.5292,'kJ/mol'), semiclassical=False),
    ],
    spinMultiplicity = 6,
    opticalIsomers = 1,
    molecularWeight = (91.0874,'amu'),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[1.03078,0.0619084,-6.49533e-05,3.60553e-08,-7.95264e-12,126891,27.0803], Tmin=(100,'K'), Tmax=(1106.58,'K')), NASAPolynomial(coeffs=[12.8376,0.0192298,-7.10116e-06,1.20191e-09,-7.85134e-14,124278,-31.0766], Tmin=(1106.58,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(1054.11,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(220.334,'J/(mol*K)'), comment="""Thermo group additivity estimation: group(O2s-(Cds-Cd)H) + group(Cs-(Cds-Cds)CtHH) + group(Cs-(Cds-Cds)HHH) + group(Cds-CdsCsCs) + group(Cds-CdsOsH) + group(Ct-CtCs) + group(Ct-CtH) + radical(AllylJ2_triplet) + radical(C=COJ) + radical(Cs_S) + radical(Acetyl)"""),
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
    label = 'C#CC=C=C[O](23376)',
    structure = SMILES('C#CC=C=C[O]'),
    E0 = (343.956,'kJ/mol'),
    modes = [
        HarmonicOscillator(frequencies=([750,770,3400,2100,2995,3025,975,1000,1300,1375,400,500,1630,1680,2175,525,540,610,2055,180],'cm^-1')),
        HinderedRotor(inertia=(1.6143,'amu*angstrom^2'), symmetry=1, barrier=(37.1158,'kJ/mol'), semiclassical=False),
    ],
    spinMultiplicity = 2,
    opticalIsomers = 1,
    molecularWeight = (79.0767,'amu'),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[1.58392,0.0406745,-9.28843e-06,-2.95441e-08,1.76387e-11,41466.8,17.6904], Tmin=(100,'K'), Tmax=(933.363,'K')), NASAPolynomial(coeffs=[17.8115,0.00285504,5.05668e-07,-1.22779e-10,4.01942e-15,37055.7,-66.882], Tmin=(933.363,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(343.956,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(203.705,'J/(mol*K)'), comment="""Thermo group additivity estimation: group(O2s-(Cds-Cd)H) + group(Cds-CdsCtH) + group(Cds-CdsOsH) + group(Cdd-CdsCds) + group(Ct-Ct(Cds-Cds)) + group(Ct-CtH) + radical(C=COJ)"""),
)

species(
    label = '[CH]C(=C=O)C[C]=[CH](24774)',
    structure = SMILES('[CH]C(=C=O)C[C]=[CH]'),
    E0 = (827.452,'kJ/mol'),
    modes = [
        HarmonicOscillator(frequencies=([3120,650,792.5,1650,2750,2850,1437.5,1250,1305,750,350,2120,512.5,787.5,350,440,435,1725,1685,370,324.724,324.967,325.005,325.326],'cm^-1')),
        HinderedRotor(inertia=(0.701109,'amu*angstrom^2'), symmetry=1, barrier=(52.5135,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.700686,'amu*angstrom^2'), symmetry=1, barrier=(52.5102,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.701946,'amu*angstrom^2'), symmetry=1, barrier=(52.5124,'kJ/mol'), semiclassical=False),
    ],
    spinMultiplicity = 5,
    opticalIsomers = 1,
    molecularWeight = (92.0954,'amu'),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[0.931706,0.0751969,-0.000115634,1.0236e-07,-3.56652e-11,99622.6,26.1391], Tmin=(100,'K'), Tmax=(850.909,'K')), NASAPolynomial(coeffs=[6.37639,0.034153,-1.60471e-05,2.9987e-09,-2.03404e-13,99255.3,4.03718], Tmin=(850.909,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(827.452,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(245.277,'J/(mol*K)'), comment="""Thermo group additivity estimation: group(Cs-Cd(CCO)HH) + group(Cs-(Cds-Cdd-O2d)HHH) + group(Cds-(Cdd-O2d)CsCs) + other(ketene_1C-C_1C-H) + group(Cds-CdsCsH) + group(Cds-CdsHH) + radical(Cds_S) + radical(AllylJ2_triplet) + radical(Cds_P)"""),
)

species(
    label = '[C]=[C]CC(=[CH])C=O(24783)',
    structure = SMILES('[C]=[C]CC(=[CH])C=O'),
    E0 = (1000.98,'kJ/mol'),
    modes = [
        HarmonicOscillator(frequencies=([1685,370,2750,2850,1437.5,1250,1305,750,350,350,440,435,1725,3120,650,792.5,1650,2782.5,750,1395,475,1775,1000,356.551],'cm^-1')),
        HinderedRotor(inertia=(0.43025,'amu*angstrom^2'), symmetry=1, barrier=(11.137,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.214523,'amu*angstrom^2'), symmetry=1, barrier=(11.1186,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.124242,'amu*angstrom^2'), symmetry=1, barrier=(11.143,'kJ/mol'), semiclassical=False),
    ],
    spinMultiplicity = 5,
    opticalIsomers = 1,
    molecularWeight = (92.0954,'amu'),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[1.17753,0.0692018,-0.000110178,1.0076e-07,-3.72108e-11,120485,26.8307], Tmin=(100,'K'), Tmax=(766.693,'K')), NASAPolynomial(coeffs=[6.77533,0.0306335,-1.64015e-05,3.28949e-09,-2.33855e-13,119902,3.10642], Tmin=(766.693,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(1000.98,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(245.277,'J/(mol*K)'), comment="""Thermo group additivity estimation: group(Cs-(Cds-Cds)(Cds-Cds)HH) + group(Cd-CdCs(CO)) + group(Cds-CdsCsH) + group(Cds-O2d(Cds-Cds)H) + group(Cds-CdsHH) + group(Cds-CdsHH) + radical(CdCdJ2_triplet) + radical(Cds_S) + radical(Cds_P)"""),
)

species(
    label = '[C]#CC=C([CH])[CH]O(26134)',
    structure = SMILES('[C]#CC=C([CH])[CH]O'),
    E0 = (867.614,'kJ/mol'),
    modes = [
        HarmonicOscillator(frequencies=([3025,407.5,1350,352.5,3010,987.5,1337.5,450,1655,3615,1277.5,1000,2175,525,350,440,435,1725,200,800,1066.67,1333.33,1600],'cm^-1')),
        HinderedRotor(inertia=(0.156089,'amu*angstrom^2'), symmetry=1, barrier=(3.5888,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.156089,'amu*angstrom^2'), symmetry=1, barrier=(3.5888,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.156089,'amu*angstrom^2'), symmetry=1, barrier=(3.5888,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.156089,'amu*angstrom^2'), symmetry=1, barrier=(3.5888,'kJ/mol'), semiclassical=False),
    ],
    spinMultiplicity = 5,
    opticalIsomers = 1,
    molecularWeight = (92.0954,'amu'),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[0.930756,0.0601245,-4.58349e-05,1.20117e-08,1.0771e-12,104467,26.671], Tmin=(100,'K'), Tmax=(986.071,'K')), NASAPolynomial(coeffs=[14.5309,0.0198633,-7.26773e-06,1.26912e-09,-8.64581e-14,101060,-42.4262], Tmin=(986.071,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(867.614,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(241.12,'J/(mol*K)'), comment="""Thermo group additivity estimation: group(O2s-CsH) + group(Cs-(Cds-Cds)OsHH) + group(Cs-(Cds-Cds)HHH) + group(Cds-CdsCsCs) + group(Cds-CdsCtH) + group(Ct-Ct(Cds-Cds)) + group(Ct-CtH) + radical(AllylJ2_triplet) + radical(Acetyl) + radical(C=CCJO)"""),
)

species(
    label = '[CH][C]1C([O])C1C#C(26135)',
    structure = SMILES('[CH][C]1C([O])C1C#C'),
    E0 = (905.13,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
    molecularWeight = (92.0954,'amu'),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[1.61765,0.0373514,1.60267e-05,-5.80893e-08,2.82226e-11,108962,22.6272], Tmin=(100,'K'), Tmax=(921.608,'K')), NASAPolynomial(coeffs=[17.1974,0.00873351,-8.74563e-07,5.6031e-11,-6.57879e-15,104434,-60.2505], Tmin=(921.608,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(905.13,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(249.434,'J/(mol*K)'), comment="""Thermo group additivity estimation: group(O2s-CsH) + group(Cs-CsCsCsH) + group(Cs-CtCsCsH) + group(Cs-CsCsOsH) + group(Cs-CsHHH) + group(Ct-CtCs) + group(Ct-CtH) + ring(Cyclopropane) + radical(CCJ2_triplet) + radical(CC(C)OJ) + radical(CCJ(C)CO)"""),
)

species(
    label = '[CH]C1[CH]OC=[C]C=1(25326)',
    structure = SMILES('[CH]C1[CH]OC=[C]C=1'),
    E0 = (590.289,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
    molecularWeight = (92.0954,'amu'),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[1.49696,0.0374521,3.27654e-05,-7.66542e-08,3.42914e-11,71101.6,18.0905], Tmin=(100,'K'), Tmax=(934.673,'K')), NASAPolynomial(coeffs=[16.9932,0.0157553,-4.02357e-06,6.61802e-10,-4.99408e-14,66255.8,-66.0495], Tmin=(934.673,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(590.289,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(253.591,'J/(mol*K)'), comment="""Thermo group additivity estimation: group(O2s-Cs(Cds-Cd)) + group(Cs-(Cds-Cds)OsHH) + group(Cs-(Cds-Cds)HHH) + group(Cds-CdsCsCs) + group(Cds-Cds(Cds-Cds)H) + group(Cds-Cds(Cds-Cds)H) + group(Cds-CdsOsH) + ring(1,3-Cyclohexadiene) + radical(C=CCJ(O)C) + radical(C=CJC=C) + radical(AllylJ2_triplet)"""),
)

species(
    label = '[CH]C(=C=O)CC#C(26136)',
    structure = SMILES('[CH]C(=C=O)CC#C'),
    E0 = (508.465,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
    molecularWeight = (92.0954,'amu'),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[0.980921,0.0703099,-9.4215e-05,7.35067e-08,-2.30587e-11,61259.1,22.7038], Tmin=(100,'K'), Tmax=(876.924,'K')), NASAPolynomial(coeffs=[8.42772,0.0293814,-1.22991e-05,2.17968e-09,-1.43727e-13,60220.7,-10.7187], Tmin=(876.924,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(508.465,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(245.277,'J/(mol*K)'), comment="""Thermo group additivity estimation: group(Cs-Cd(CCO)HH) + group(Cs-(Cds-Cdd-O2d)HHH) + group(Cds-(Cdd-O2d)CsCs) + other(ketene_1C-C_1C-H) + group(Ct-CtCs) + group(Ct-CtH) + radical(AllylJ2_triplet)"""),
)

species(
    label = '[CH]C1=COC1C#C(25333)',
    structure = SMILES('[CH]C1=COC1C#C'),
    E0 = (560.328,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
    molecularWeight = (92.0954,'amu'),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[1.07065,0.0413608,3.93092e-05,-9.92228e-08,4.65717e-11,67518.8,20.3307], Tmin=(100,'K'), Tmax=(915.724,'K')), NASAPolynomial(coeffs=[23.5135,0.00356878,2.53598e-06,-6.11076e-10,3.70404e-14,60882.7,-99.7586], Tmin=(915.724,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(560.328,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(249.434,'J/(mol*K)'), comment="""Thermo group additivity estimation: group(O2s-Cs(Cds-Cd)) + group(Cs-(Cds-Cds)(Cds-Cds)OsH) + group(Cs-(Cds-Cds)HHH) + group(Cds-CdsCsCs) + group(Cds-CdsOsH) + group(Ct-CtCs) + group(Ct-CtH) + ring(Cyclobutene) + radical(AllylJ2_triplet)"""),
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
    label = '[CH]C([CH])=C[O](21383)',
    structure = SMILES('[CH]C([CH])=C[O]'),
    E0 = (641.108,'kJ/mol'),
    modes = [
        HarmonicOscillator(frequencies=([350,440,435,1725,3010,987.5,1337.5,450,1655,180,180,180,180,180,180,180],'cm^-1')),
        HinderedRotor(inertia=(2.17338,'amu*angstrom^2'), symmetry=1, barrier=(49.9703,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(2.17106,'amu*angstrom^2'), symmetry=1, barrier=(49.917,'kJ/mol'), semiclassical=False),
    ],
    spinMultiplicity = 6,
    opticalIsomers = 1,
    molecularWeight = (67.066,'amu'),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[1.87426,0.0365132,4.33688e-06,-3.49083e-08,1.70988e-11,77193.3,19.0396], Tmin=(100,'K'), Tmax=(940.923,'K')), NASAPolynomial(coeffs=[12.3389,0.0174493,-5.80035e-06,9.89601e-10,-6.87008e-14,74098.6,-36.7892], Tmin=(940.923,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(641.108,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(174.604,'J/(mol*K)'), comment="""Thermo group additivity estimation: group(O2s-(Cds-Cd)H) + group(Cs-(Cds-Cds)HHH) + group(Cs-(Cds-Cds)HHH) + group(Cds-CdsCsCs) + group(Cds-CdsOsH) + radical(C=COJ) + radical(AllylJ2_triplet) + radical(AllylJ2_triplet)"""),
)

species(
    label = '[CH]C(=C=C=[CH])C=O(26137)',
    structure = SMILES('[CH]C(=C=C=[CH])C=O'),
    E0 = (679.537,'kJ/mol'),
    modes = [
        HarmonicOscillator(frequencies=([3120,650,792.5,1650,2782.5,750,1395,475,1775,1000,350,440,435,1725,540,563.333,586.667,610,1970,2140,180,180],'cm^-1')),
        HinderedRotor(inertia=(2.02811,'amu*angstrom^2'), symmetry=1, barrier=(46.6302,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(2.02746,'amu*angstrom^2'), symmetry=1, barrier=(46.6153,'kJ/mol'), semiclassical=False),
    ],
    spinMultiplicity = 4,
    opticalIsomers = 1,
    molecularWeight = (91.0874,'amu'),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[1.1304,0.0651432,-7.17438e-05,4.22118e-08,-1.01165e-11,81831.1,22.3256], Tmin=(100,'K'), Tmax=(1004.18,'K')), NASAPolynomial(coeffs=[11.2609,0.0247905,-1.14675e-05,2.19551e-09,-1.54161e-13,79796.5,-26.5907], Tmin=(1004.18,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(679.537,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(224.491,'J/(mol*K)'), comment="""Thermo group additivity estimation: group(Cs-(Cds-Cds)HHH) + group(Cds-Cds(Cds-Cds)Cs) + group(Cds-O2d(Cds-Cds)H) + group(Cds-CdsHH) + group(Cdd-CdsCds) + group(Cdd-CdsCds) + radical(AllylJ2_triplet) + radical(C=C=CJ)"""),
)

species(
    label = '[CH]C([C]=C=[CH])C=O(26138)',
    structure = SMILES('[CH]C([C]=C=[CH])C=O'),
    E0 = (851.116,'kJ/mol'),
    modes = [
        HarmonicOscillator(frequencies=([3120,650,792.5,1650,2782.5,750,1395,475,1775,1000,540,610,2055,1380,1390,370,380,2900,435,1685,370,251.321,251.973,2123.67],'cm^-1')),
        HinderedRotor(inertia=(0.147952,'amu*angstrom^2'), symmetry=1, barrier=(6.73769,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.148615,'amu*angstrom^2'), symmetry=1, barrier=(6.73762,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(1.50478,'amu*angstrom^2'), symmetry=1, barrier=(65.7424,'kJ/mol'), semiclassical=False),
    ],
    spinMultiplicity = 5,
    opticalIsomers = 1,
    molecularWeight = (92.0954,'amu'),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[0.987107,0.0724715,-0.00011408,9.67498e-08,-3.18329e-11,102468,27.3271], Tmin=(100,'K'), Tmax=(876.951,'K')), NASAPolynomial(coeffs=[8.63199,0.0249406,-1.11238e-05,2.0191e-09,-1.34005e-13,101614,-5.7756], Tmin=(876.951,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(851.116,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(245.277,'J/(mol*K)'), comment="""Thermo group additivity estimation: group(Cs-(Cds-O2d)CsCsH) + group(Cs-CsHHH) + group(Cds-CdsCsH) + group(Cds-OdCsH) + group(Cds-CdsHH) + group(Cdd-CdsCds) + radical(C=C=CJ) + radical(Cds_S) + radical(CCJ2_triplet)"""),
)

species(
    label = '[CH]C([C]=O)C=C=[CH](26139)',
    structure = SMILES('[CH]C([C]=O)C=C=[CH]'),
    E0 = (771.962,'kJ/mol'),
    modes = [
        HarmonicOscillator(frequencies=([1855,455,950,3120,650,792.5,1650,540,610,2055,1380,1390,370,380,2900,435,3010,987.5,1337.5,450,1655,339.987,340.811,2159.78],'cm^-1')),
        HinderedRotor(inertia=(0.120086,'amu*angstrom^2'), symmetry=1, barrier=(9.82298,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.119253,'amu*angstrom^2'), symmetry=1, barrier=(9.81712,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.83109,'amu*angstrom^2'), symmetry=1, barrier=(68.6211,'kJ/mol'), semiclassical=False),
    ],
    spinMultiplicity = 5,
    opticalIsomers = 1,
    molecularWeight = (92.0954,'amu'),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[1.20536,0.0648203,-8.40148e-05,5.70927e-08,-1.46308e-11,92943.3,27.2397], Tmin=(100,'K'), Tmax=(728.551,'K')), NASAPolynomial(coeffs=[10.3769,0.021439,-9.05572e-06,1.63904e-09,-1.10384e-13,91421.8,-15.3731], Tmin=(728.551,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(771.962,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(245.277,'J/(mol*K)'), comment="""Thermo group additivity estimation: group(Cs-(Cds-O2d)CsCsH) + group(Cs-CsHHH) + group(Cds-CdsCsH) + group(Cds-OdCsH) + group(Cds-CdsHH) + group(Cdd-CdsCds) + radical(C=C=CJ) + radical(CC(C)CJ=O) + radical(CCJ2_triplet)"""),
)

species(
    label = '[CH]C1=CC1([CH])C=O(26140)',
    structure = SMILES('[CH]C1=CC1([CH])C=O'),
    E0 = (910.942,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
    molecularWeight = (92.0954,'amu'),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[1.17309,0.0491906,-8.7454e-06,-2.69024e-08,1.41706e-11,109674,24.6569], Tmin=(100,'K'), Tmax=(1004.65,'K')), NASAPolynomial(coeffs=[15.9772,0.0194684,-7.99618e-06,1.55102e-09,-1.13988e-13,105225,-54.1723], Tmin=(1004.65,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(910.942,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(245.277,'J/(mol*K)'), comment="""Thermo group additivity estimation: group(Cs-(Cds-Cds)(Cds-Cds)(Cds-Cds)Cs) + group(Cs-CsHHH) + group(Cs-(Cds-Cds)HHH) + group(Cds-CdsCsCs) + group(Cds-CdsCsH) + group(Cds-OdCsH) + ring(Cyclopropene) + radical(AllylJ2_triplet) + radical(CCJ2_triplet)"""),
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
    label = '[CH][C]=C[C]=[CH](21426)',
    structure = SMILES('[CH][C]=C[C]=[CH]'),
    E0 = (1112.51,'kJ/mol'),
    modes = [
        HarmonicOscillator(frequencies=([1670,1700,300,440,3010,987.5,1337.5,450,1655,3120,650,792.5,1650,180,180,180],'cm^-1')),
        HinderedRotor(inertia=(2.13683,'amu*angstrom^2'), symmetry=1, barrier=(49.1299,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(2.1316,'amu*angstrom^2'), symmetry=1, barrier=(49.0098,'kJ/mol'), semiclassical=False),
    ],
    spinMultiplicity = 6,
    opticalIsomers = 1,
    molecularWeight = (63.0773,'amu'),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[1.95039,0.0458062,-4.82344e-05,2.71146e-08,-5.46777e-12,133878,18.0544], Tmin=(100,'K'), Tmax=(785.352,'K')), NASAPolynomial(coeffs=[8.34208,0.0187805,-7.17602e-06,1.22516e-09,-7.99173e-14,132703,-12.3229], Tmin=(785.352,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(1112.51,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(174.604,'J/(mol*K)'), comment="""Thermo group additivity estimation: group(Cs-(Cds-Cds)HHH) + group(Cds-CdsCsH) + group(Cds-Cds(Cds-Cds)H) + group(Cds-Cds(Cds-Cds)H) + group(Cds-CdsHH) + radical(Cds_S) + radical(Cds_P) + radical(C=CJC=C) + radical(AllylJ2_triplet)"""),
)

species(
    label = '[CH]=C1[CH]C(=[CH])C1[O](26113)',
    structure = SMILES('[CH]=C1[CH]C(=[CH])C1[O]'),
    E0 = (844.581,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
    molecularWeight = (92.0954,'amu'),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[1.98009,0.0219216,6.6602e-05,-1.12475e-07,4.73297e-11,101673,23.635], Tmin=(100,'K'), Tmax=(937.982,'K')), NASAPolynomial(coeffs=[19.0915,0.00626079,-3.07126e-09,3.84935e-12,-1.08299e-14,95941.9,-71.2613], Tmin=(937.982,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(844.581,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(257.749,'J/(mol*K)'), comment="""Thermo group additivity estimation: group(O2s-CsH) + group(Cs-(Cds-Cds)(Cds-Cds)OsH) + group(Cs-(Cds-Cds)(Cds-Cds)HH) + group(Cds-CdsCsCs) + group(Cds-CdsCsCs) + group(Cds-CdsHH) + group(Cds-CdsHH) + ring(Cyclobutane) + radical(Cds_P) + radical(C=CCJC=C) + radical(CC(C)OJ) + radical(Cds_P)"""),
)

species(
    label = '[CH]C=C=C([CH])C=O(26141)',
    structure = SMILES('[CH]C=C=C([CH])C=O'),
    E0 = (719.141,'kJ/mol'),
    modes = [
        HarmonicOscillator(frequencies=([2782.5,750,1395,475,1775,1000,540,610,2055,350,440,435,1725,3010,987.5,1337.5,450,1655,180,180,180,180,180,180],'cm^-1')),
        HinderedRotor(inertia=(2.13674,'amu*angstrom^2'), symmetry=1, barrier=(49.1278,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(2.13605,'amu*angstrom^2'), symmetry=1, barrier=(49.112,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(2.13517,'amu*angstrom^2'), symmetry=1, barrier=(49.0917,'kJ/mol'), semiclassical=False),
    ],
    spinMultiplicity = 5,
    opticalIsomers = 1,
    molecularWeight = (92.0954,'amu'),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[1.27014,0.0623114,-4.74781e-05,1.87473e-08,-3.154e-12,86588.6,23.2403], Tmin=(100,'K'), Tmax=(1330.82,'K')), NASAPolynomial(coeffs=[10.4805,0.0346282,-1.62755e-05,3.11648e-09,-2.17682e-13,84137.2,-23.8264], Tmin=(1330.82,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(719.141,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(245.277,'J/(mol*K)'), comment="""Thermo group additivity estimation: group(Cs-(Cds-Cds)HHH) + group(Cs-(Cds-Cds)HHH) + group(Cds-Cds(Cds-Cds)Cs) + group(Cds-CdsCsH) + group(Cds-O2d(Cds-Cds)H) + group(Cdd-CdsCds) + radical(AllylJ2_triplet) + radical(AllylJ2_triplet)"""),
)

species(
    label = '[CH]C(=C=O)[CH]C=[CH](26142)',
    structure = SMILES('[CH]C(=C=O)[CH]C=[CH]'),
    E0 = (690.741,'kJ/mol'),
    modes = [
        HarmonicOscillator(frequencies=([3120,650,792.5,1650,3025,407.5,1350,352.5,2120,512.5,787.5,350,440,435,1725,3010,987.5,1337.5,450,1655,325.529,325.529,325.53,325.53],'cm^-1')),
        HinderedRotor(inertia=(0.674703,'amu*angstrom^2'), symmetry=1, barrier=(50.7366,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.674712,'amu*angstrom^2'), symmetry=1, barrier=(50.7366,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.674702,'amu*angstrom^2'), symmetry=1, barrier=(50.7366,'kJ/mol'), semiclassical=False),
    ],
    spinMultiplicity = 5,
    opticalIsomers = 1,
    molecularWeight = (92.0954,'amu'),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[1.25552,0.0606918,-5.65402e-05,2.88878e-08,-6.1378e-12,83175.6,23.026], Tmin=(100,'K'), Tmax=(1113.48,'K')), NASAPolynomial(coeffs=[10.4525,0.0276531,-1.20329e-05,2.24022e-09,-1.54846e-13,81127.4,-22.3327], Tmin=(1113.48,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(690.741,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(245.277,'J/(mol*K)'), comment="""Thermo group additivity estimation: group(Cs-Cd(CCO)HH) + group(Cs-(Cds-Cdd-O2d)HHH) + group(Cds-(Cdd-O2d)CsCs) + other(ketene_1C-C_1C-H) + group(Cds-CdsCsH) + group(Cds-CdsHH) + radical(C=CCJC=C=O) + radical(AllylJ2_triplet) + radical(Cds_P)"""),
)

species(
    label = '[CH][C]=CC([CH])[C]=O(26143)',
    structure = SMILES('[CH][C]=CC([CH])[C]=O'),
    E0 = (1049.41,'kJ/mol'),
    modes = [
        HarmonicOscillator(frequencies=([1855,455,950,1685,370,1380,1390,370,380,2900,435,3010,987.5,1337.5,450,1655,200,800,960,1120,1280,1440,1600],'cm^-1')),
        HinderedRotor(inertia=(0.156089,'amu*angstrom^2'), symmetry=1, barrier=(3.5888,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.156089,'amu*angstrom^2'), symmetry=1, barrier=(3.5888,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.156089,'amu*angstrom^2'), symmetry=1, barrier=(3.5888,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.156089,'amu*angstrom^2'), symmetry=1, barrier=(3.5888,'kJ/mol'), semiclassical=False),
    ],
    spinMultiplicity = 7,
    opticalIsomers = 1,
    molecularWeight = (92.0954,'amu'),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[1.08385,0.0693376,-9.48732e-05,7.83037e-08,-2.65717e-11,126315,29.5664], Tmin=(100,'K'), Tmax=(806.621,'K')), NASAPolynomial(coeffs=[7.16745,0.0328704,-1.53448e-05,2.89279e-09,-1.9875e-13,125538,2.79409], Tmin=(806.621,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(1049.41,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(241.12,'J/(mol*K)'), comment="""Thermo group additivity estimation: group(Cs-(Cds-O2d)CsCsH) + group(Cs-CsHHH) + group(Cs-(Cds-Cds)HHH) + group(Cds-CdsCsH) + group(Cds-CdsCsH) + group(Cds-OdCsH) + radical(CCJ2_triplet) + radical(Cds_S) + radical(CC(C)CJ=O) + radical(AllylJ2_triplet)"""),
)

species(
    label = '[CH]C(=[C][C]=[CH])C[O](26144)',
    structure = SMILES('[CH]C(=[C][C]=[CH])C[O]'),
    E0 = (1107.18,'kJ/mol'),
    modes = [
        HarmonicOscillator(frequencies=([3120,650,792.5,1650,2750,2850,1437.5,1250,1305,750,350,350,440,435,1725,1670,1700,300,440,285.793,285.795,285.796,285.797,285.797],'cm^-1')),
        HinderedRotor(inertia=(0.909585,'amu*angstrom^2'), symmetry=1, barrier=(52.7213,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.909605,'amu*angstrom^2'), symmetry=1, barrier=(52.7213,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.909599,'amu*angstrom^2'), symmetry=1, barrier=(52.7213,'kJ/mol'), semiclassical=False),
    ],
    spinMultiplicity = 7,
    opticalIsomers = 1,
    molecularWeight = (92.0954,'amu'),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[0.783506,0.0807398,-0.000130789,1.17974e-07,-4.05263e-11,133270,27.4244], Tmin=(100,'K'), Tmax=(898.918,'K')), NASAPolynomial(coeffs=[5.40291,0.0360197,-1.58424e-05,2.82079e-09,-1.83996e-13,133415,11.0613], Tmin=(898.918,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(1107.18,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(245.277,'J/(mol*K)'), comment="""Thermo group additivity estimation: group(O2s-CsH) + group(Cs-(Cds-Cds)OsHH) + group(Cs-(Cds-Cds)HHH) + group(Cds-CdsCsCs) + group(Cds-Cds(Cds-Cds)H) + group(Cds-Cds(Cds-Cds)H) + group(Cds-CdsHH) + radical(C=CJC=C) + radical(Cds_P) + radical(AllylJ2_triplet) + radical(C=CJC=C) + radical(CCOJ)"""),
)

species(
    label = '[CH][C]=[C]C([CH])=CO(26145)',
    structure = SMILES('[CH][C]=[C]C([CH])=CO'),
    E0 = (989.718,'kJ/mol'),
    modes = [
        HarmonicOscillator(frequencies=([1670,1700,300,440,3010,987.5,1337.5,450,1655,3615,1277.5,1000,350,440,435,1725,200,800,960,1120,1280,1440,1600],'cm^-1')),
        HinderedRotor(inertia=(0.156089,'amu*angstrom^2'), symmetry=1, barrier=(3.5888,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.156089,'amu*angstrom^2'), symmetry=1, barrier=(3.5888,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.156089,'amu*angstrom^2'), symmetry=1, barrier=(3.5888,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.156089,'amu*angstrom^2'), symmetry=1, barrier=(3.5888,'kJ/mol'), semiclassical=False),
    ],
    spinMultiplicity = 7,
    opticalIsomers = 1,
    molecularWeight = (92.0954,'amu'),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[-0.437821,0.0850118,-9.18142e-05,4.97331e-08,-1.02647e-11,119206,28.472], Tmin=(100,'K'), Tmax=(1324.59,'K')), NASAPolynomial(coeffs=[19.652,0.0156889,-3.50937e-06,3.55937e-10,-1.42677e-14,114643,-71.2307], Tmin=(1324.59,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(989.718,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(241.12,'J/(mol*K)'), comment="""Thermo group additivity estimation: group(O2s-(Cds-Cd)H) + group(Cs-(Cds-Cds)HHH) + group(Cs-(Cds-Cds)HHH) + group(Cds-Cds(Cds-Cds)Cs) + group(Cds-CdsCsH) + group(Cds-Cds(Cds-Cds)H) + group(Cds-CdsOsH) + radical(AllylJ2_triplet) + radical(C=CJC=C) + radical(Cds_S) + radical(AllylJ2_triplet)"""),
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
    label = '[CH][CH]C=C=[CH](19252)',
    structure = SMILES('[CH][CH]C=C=[CH]'),
    E0 = (867.861,'kJ/mol'),
    modes = [
        HarmonicOscillator(frequencies=([3120,650,792.5,1650,3025,407.5,1350,352.5,3010,987.5,1337.5,450,1655,540,610,2055,524.185,563.719,564.949],'cm^-1')),
        HinderedRotor(inertia=(0.0200302,'amu*angstrom^2'), symmetry=1, barrier=(79.5804,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.354495,'amu*angstrom^2'), symmetry=1, barrier=(79.5766,'kJ/mol'), semiclassical=False),
    ],
    spinMultiplicity = 5,
    opticalIsomers = 1,
    molecularWeight = (64.0853,'amu'),
    collisionModel = TransportData(shapeIndex=2, epsilon=(3149.55,'J/mol'), sigma=(5.50528,'angstroms'), dipoleMoment=(0,'C*m'), polarizability=(0,'angstroms^3'), rotrelaxcollnum=0, comment="""Epsilon & sigma estimated with Tc=491.95 K, Pc=42.83 bar (from Joback method)"""),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[2.19239,0.0344402,-1.76598e-05,-4.40435e-09,5.18614e-12,104450,17.9495], Tmin=(100,'K'), Tmax=(942.255,'K')), NASAPolynomial(coeffs=[10.2948,0.0130216,-4.22214e-06,7.05005e-10,-4.76148e-14,102347,-23.7152], Tmin=(942.255,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(867.861,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(199.547,'J/(mol*K)'), comment="""Thermo group additivity estimation: group(Cs-(Cds-Cds)CsHH) + group(Cs-CsHHH) + group(Cds-CdsCsH) + group(Cds-CdsHH) + group(Cdd-CdsCds) + radical(C=C=CJ) + radical(CCJ2_triplet) + radical(Allyl_S)"""),
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
    label = '[CH]C(=C=C=[CH])C[O](26146)',
    structure = SMILES('[CH]C(=C=C=[CH])C[O]'),
    E0 = (845.957,'kJ/mol'),
    modes = [
        HarmonicOscillator(frequencies=([3120,650,792.5,1650,2750,2850,1437.5,1250,1305,750,350,350,440,435,1725,540,563.333,586.667,610,1970,2140,244.167,244.171,244.176,244.177],'cm^-1')),
        HinderedRotor(inertia=(1.23118,'amu*angstrom^2'), symmetry=1, barrier=(52.0886,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(1.23119,'amu*angstrom^2'), symmetry=1, barrier=(52.0887,'kJ/mol'), semiclassical=False),
    ],
    spinMultiplicity = 5,
    opticalIsomers = 1,
    molecularWeight = (92.0954,'amu'),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[1.11855,0.0709603,-0.000104888,9.32226e-08,-3.26878e-11,101842,26.3355], Tmin=(100,'K'), Tmax=(862.959,'K')), NASAPolynomial(coeffs=[5.02461,0.0365064,-1.65834e-05,3.05066e-09,-2.0505e-13,101776,11.5937], Tmin=(862.959,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(845.957,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(249.434,'J/(mol*K)'), comment="""Thermo group additivity estimation: group(O2s-CsH) + group(Cs-(Cds-Cds)OsHH) + group(Cs-(Cds-Cds)HHH) + group(Cds-CdsCsCs) + group(Cds-CdsHH) + group(Cdd-CdsCds) + group(Cdd-CdsCds) + radical(CCOJ) + radical(AllylJ2_triplet) + radical(C=C=CJ)"""),
)

species(
    label = '[C]#CC=C([CH])C[O](26147)',
    structure = SMILES('[C]#CC=C([CH])C[O]'),
    E0 = (976.024,'kJ/mol'),
    modes = [
        HarmonicOscillator(frequencies=([2750,2850,1437.5,1250,1305,750,350,350,440,435,1725,2175,525,3010,987.5,1337.5,450,1655,445.387,445.39,445.391,445.391,445.391,445.392],'cm^-1')),
        HinderedRotor(inertia=(0.375026,'amu*angstrom^2'), symmetry=1, barrier=(52.7938,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.375018,'amu*angstrom^2'), symmetry=1, barrier=(52.7937,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.375039,'amu*angstrom^2'), symmetry=1, barrier=(52.7938,'kJ/mol'), semiclassical=False),
    ],
    spinMultiplicity = 5,
    opticalIsomers = 1,
    molecularWeight = (92.0954,'amu'),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[1.26019,0.0663637,-9.09332e-05,7.86536e-08,-2.77114e-11,117481,25.7532], Tmin=(100,'K'), Tmax=(834.493,'K')), NASAPolynomial(coeffs=[5.14739,0.0364091,-1.67384e-05,3.1219e-09,-2.12647e-13,117227,10.0651], Tmin=(834.493,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(976.024,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(245.277,'J/(mol*K)'), comment="""Thermo group additivity estimation: group(O2s-CsH) + group(Cs-(Cds-Cds)OsHH) + group(Cs-(Cds-Cds)HHH) + group(Cds-CdsCsCs) + group(Cds-CdsCtH) + group(Ct-Ct(Cds-Cds)) + group(Ct-CtH) + radical(CCOJ) + radical(AllylJ2_triplet) + radical(Acetyl)"""),
)

species(
    label = '[C]#CC[C]([CH])[CH][O](26148)',
    structure = SMILES('[C]#CC[C]([CH])[CH][O]'),
    E0 = (1277.09,'kJ/mol'),
    modes = [
        HarmonicOscillator(frequencies=([3025,407.5,1350,352.5,2750,2850,1437.5,1250,1305,750,350,2175,525,360,370,350,200,800,960,1120,1280,1440,1600],'cm^-1')),
        HinderedRotor(inertia=(0.156089,'amu*angstrom^2'), symmetry=1, barrier=(3.5888,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.156089,'amu*angstrom^2'), symmetry=1, barrier=(3.5888,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.156089,'amu*angstrom^2'), symmetry=1, barrier=(3.5888,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.156089,'amu*angstrom^2'), symmetry=1, barrier=(3.5888,'kJ/mol'), semiclassical=False),
    ],
    spinMultiplicity = 7,
    opticalIsomers = 1,
    molecularWeight = (92.0954,'amu'),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[0.954317,0.0765328,-0.000130484,1.17322e-07,-3.96776e-11,153699,28.1169], Tmin=(100,'K'), Tmax=(907.701,'K')), NASAPolynomial(coeffs=[6.54472,0.0284874,-1.24024e-05,2.18299e-09,-1.40495e-13,153648,7.00015], Tmin=(907.701,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(1277.09,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(241.12,'J/(mol*K)'), comment="""Thermo group additivity estimation: group(O2s-CsH) + group(Cs-CsCsCsH) + group(Cs-CtCsHH) + group(Cs-CsOsHH) + group(Cs-CsHHH) + group(Ct-CtCs) + group(Ct-CtH) + radical(Acetyl) + radical(CCJ2_triplet) + radical(CCOJ) + radical(CCsJOH) + radical(CCJ(C)CO)"""),
)

species(
    label = '[C]#C[CH]C([CH])[CH][O](26149)',
    structure = SMILES('[C]#C[CH]C([CH])[CH][O]'),
    E0 = (1270.72,'kJ/mol'),
    modes = [
        HarmonicOscillator(frequencies=([3000,3050,390,425,1340,1360,335,370,1380,1390,370,380,2900,435,2175,525,200,800,960,1120,1280,1440,1600],'cm^-1')),
        HinderedRotor(inertia=(0.156089,'amu*angstrom^2'), symmetry=1, barrier=(3.5888,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.156089,'amu*angstrom^2'), symmetry=1, barrier=(3.5888,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.156089,'amu*angstrom^2'), symmetry=1, barrier=(3.5888,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.156089,'amu*angstrom^2'), symmetry=1, barrier=(3.5888,'kJ/mol'), semiclassical=False),
    ],
    spinMultiplicity = 7,
    opticalIsomers = 1,
    molecularWeight = (92.0954,'amu'),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[0.688243,0.0839468,-0.000150093,1.34528e-07,-4.44738e-11,152942,28.7147], Tmin=(100,'K'), Tmax=(937.153,'K')), NASAPolynomial(coeffs=[7.31464,0.0262158,-1.05557e-05,1.73427e-09,-1.04515e-13,152993,4.07556], Tmin=(937.153,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(1270.72,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(241.12,'J/(mol*K)'), comment="""Thermo group additivity estimation: group(O2s-CsH) + group(Cs-CsCsCsH) + group(Cs-CtCsHH) + group(Cs-CsOsHH) + group(Cs-CsHHH) + group(Ct-CtCs) + group(Ct-CtH) + radical(CCsJOH) + radical(CCJ2_triplet) + radical(CCOJ) + radical(Acetyl) + radical(Sec_Propargyl)"""),
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
    label = '[CH]=[C]C=C=[CH](21430)',
    structure = SMILES('[CH]=[C]C=C=[CH]'),
    E0 = (835.068,'kJ/mol'),
    modes = [
        HarmonicOscillator(frequencies=([3115,3125,620,680,785,800,1600,1700,1685,370,540,610,2055,197.988,198.111,199.477,200.18],'cm^-1')),
        HinderedRotor(inertia=(0.0653758,'amu*angstrom^2'), symmetry=1, barrier=(71.8386,'kJ/mol'), semiclassical=False),
    ],
    spinMultiplicity = 4,
    opticalIsomers = 1,
    molecularWeight = (63.0773,'amu'),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[1.98287,0.0425038,-4.25398e-05,1.42591e-08,1.96062e-12,100510,16.0374], Tmin=(100,'K'), Tmax=(798.852,'K')), NASAPolynomial(coeffs=[11.6648,0.0071383,-7.57417e-07,-6.06055e-11,1.11936e-14,98544.8,-31.1168], Tmin=(798.852,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(835.068,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(178.761,'J/(mol*K)'), comment="""Thermo group additivity estimation: group(Cds-Cds(Cds-Cds)H) + group(Cds-Cds(Cds-Cds)H) + group(Cds-CdsHH) + group(Cds-CdsHH) + group(Cdd-CdsCds) + radical(C=C=CJ) + radical(C=CJC=C) + radical(Cds_P)"""),
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
    label = '[CH]=C=CC#C(26151)',
    structure = SMILES('[CH]=C=CC#C'),
    E0 = (565.763,'kJ/mol'),
    modes = [
        HarmonicOscillator(frequencies=([3120,650,792.5,1650,3010,987.5,1337.5,450,1655,2175,525,540,610,2055,369.038,371.106,372.113],'cm^-1')),
        HinderedRotor(inertia=(0.523585,'amu*angstrom^2'), symmetry=1, barrier=(51.2195,'kJ/mol'), semiclassical=False),
    ],
    spinMultiplicity = 2,
    opticalIsomers = 1,
    molecularWeight = (63.0773,'amu'),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[2.22155,0.032884,-1.60645e-05,-8.96827e-09,7.85889e-12,68115.4,14.066], Tmin=(100,'K'), Tmax=(916.339,'K')), NASAPolynomial(coeffs=[11.802,0.00757715,-1.67028e-06,2.25937e-10,-1.50201e-14,65666.3,-35.0998], Tmin=(916.339,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(565.763,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(178.761,'J/(mol*K)'), comment="""Thermo group additivity estimation: group(Cds-CdsCtH) + group(Cds-CdsHH) + group(Cdd-CdsCds) + group(Ct-Ct(Cds-Cds)) + group(Ct-CtH) + radical(C=C=CJ)"""),
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
    label = '[CH]=C=[C]C([CH2])=C[O](26152)',
    structure = SMILES('[CH]=C=[C]C([CH2])=C[O]'),
    E0 = (673.397,'kJ/mol'),
    modes = [
        HarmonicOscillator(frequencies=([3120,650,792.5,1650,540,610,2055,350,440,435,1725,3000,3100,440,815,1455,1000,1685,370,3010,987.5,1337.5,450,1655,180],'cm^-1')),
        HinderedRotor(inertia=(1.5466,'amu*angstrom^2'), symmetry=1, barrier=(35.5594,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(1.56677,'amu*angstrom^2'), symmetry=1, barrier=(36.0232,'kJ/mol'), semiclassical=False),
    ],
    spinMultiplicity = 5,
    opticalIsomers = 1,
    molecularWeight = (92.0954,'amu'),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[-0.262936,0.0740095,-8.2781e-05,4.30509e-08,-8.19208e-12,81161.4,27.4083], Tmin=(100,'K'), Tmax=(1520.67,'K')), NASAPolynomial(coeffs=[21.0776,0.00293873,2.057e-06,-6.01383e-10,4.62868e-14,76398,-78.8138], Tmin=(1520.67,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(673.397,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(249.434,'J/(mol*K)'), comment="""Thermo group additivity estimation: group(O2s-(Cds-Cd)H) + group(Cs-(Cds-Cds)HHH) + group(Cds-Cds(Cds-Cds)Cs) + group(Cds-Cds(Cds-Cds)H) + group(Cds-CdsOsH) + group(Cds-CdsHH) + group(Cdd-CdsCds) + radical(Allyl_P) + radical(Cds_S) + radical(C=C=CJ) + radical(C=COJ)"""),
)

species(
    label = '[CH]=C1[CH]C([CH]1)=C[O](26092)',
    structure = SMILES('[CH]=C1[CH]C([CH]1)=C[O]'),
    E0 = (579.352,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
    molecularWeight = (92.0954,'amu'),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[2.27375,0.00189327,0.000148125,-2.15217e-07,8.90628e-11,69775.8,20.2337], Tmin=(100,'K'), Tmax=(914.304,'K')), NASAPolynomial(coeffs=[25.6259,-0.00674233,8.85079e-06,-1.78336e-09,1.11016e-13,61596.3,-111.713], Tmin=(914.304,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(579.352,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(257.749,'J/(mol*K)'), comment="""Thermo group additivity estimation: group(O2s-(Cds-Cd)H) + group(Cs-(Cds-Cds)(Cds-Cds)HH) + group(Cs-(Cds-Cds)(Cds-Cds)HH) + group(Cds-CdsCsCs) + group(Cds-CdsCsCs) + group(Cds-CdsOsH) + group(Cds-CdsHH) + ring(Cyclobutane) + radical(C=CCJC=C) + radical(C=CCJC=C) + radical(Cds_P) + radical(C=COJ)"""),
)

species(
    label = '[CH][C]=CC([O])C#C(22636)',
    structure = SMILES('[CH][C]=CC([O])C#C'),
    E0 = (879.728,'kJ/mol'),
    modes = [
        HarmonicOscillator(frequencies=([1685,370,750,770,3400,2100,1380,1390,370,380,2900,435,2175,525,3010,987.5,1337.5,450,1655,275.103,275.163,275.642,275.807,276.926],'cm^-1')),
        HinderedRotor(inertia=(0.928761,'amu*angstrom^2'), symmetry=1, barrier=(50.547,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.942385,'amu*angstrom^2'), symmetry=1, barrier=(50.5704,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.921895,'amu*angstrom^2'), symmetry=1, barrier=(50.5884,'kJ/mol'), semiclassical=False),
    ],
    spinMultiplicity = 5,
    opticalIsomers = 1,
    molecularWeight = (92.0954,'amu'),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[1.15544,0.062126,-6.17755e-05,3.39125e-08,-7.6334e-12,105910,26.9269], Tmin=(100,'K'), Tmax=(1066.04,'K')), NASAPolynomial(coeffs=[10.8452,0.0257676,-1.06159e-05,1.91858e-09,-1.30324e-13,103844,-20.4402], Tmin=(1066.04,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(879.728,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(245.277,'J/(mol*K)'), comment="""Thermo group additivity estimation: group(O2s-CsH) + group(Cs-(Cds-Cds)(Cds-Cds)OsH) + group(Cs-(Cds-Cds)HHH) + group(Cds-CdsCsH) + group(Cds-CdsCsH) + group(Ct-CtCs) + group(Ct-CtH) + radical(CC(C)OJ) + radical(Cds_S) + radical(AllylJ2_triplet)"""),
)

species(
    label = 'C#CC=C1[CH]C1[O](26153)',
    structure = SMILES('C#CC=C1[CH]C1[O]'),
    E0 = (589.973,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
    molecularWeight = (92.0954,'amu'),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[1.14457,0.0460817,2.23996e-06,-4.73146e-08,2.41395e-11,71075.4,19.5358], Tmin=(100,'K'), Tmax=(956.822,'K')), NASAPolynomial(coeffs=[20.1343,0.00771977,-1.93428e-06,4.04524e-10,-3.68729e-14,65563.5,-81.0539], Tmin=(956.822,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(589.973,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(253.591,'J/(mol*K)'), comment="""Thermo group additivity estimation: group(O2s-CsH) + group(Cs-(Cds-Cds)CsOsH) + group(Cs-(Cds-Cds)CsHH) + group(Cds-CdsCsCs) + group(Cds-CdsCtH) + group(Ct-Ct(Cds-Cds)) + group(Ct-CtH) + ring(Methylene_cyclopropane) + radical(CC(C)OJ) + radical(C=CCJCO)"""),
)

species(
    label = '[C][C](C=O)C=C=[CH](26154)',
    structure = SMILES('[C][C](C=O)C=C=[CH]'),
    E0 = (954.237,'kJ/mol'),
    modes = [
        HarmonicOscillator(frequencies=([3120,650,792.5,1650,360,370,350,540,610,2055,3010,987.5,1337.5,450,1655,2782.5,750,1395,475,1775,1000,180],'cm^-1')),
        HinderedRotor(inertia=(2.06823,'amu*angstrom^2'), symmetry=1, barrier=(47.5527,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(2.06569,'amu*angstrom^2'), symmetry=1, barrier=(47.4943,'kJ/mol'), semiclassical=False),
    ],
    spinMultiplicity = 6,
    opticalIsomers = 1,
    molecularWeight = (91.0874,'amu'),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[1.40464,0.05139,-4.52319e-05,1.73516e-08,-1.62979e-12,114867,22.2615], Tmin=(100,'K'), Tmax=(1003.91,'K')), NASAPolynomial(coeffs=[13.2642,0.0143636,-5.18961e-06,9.08486e-10,-6.20855e-14,111970,-37.5672], Tmin=(1003.91,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(954.237,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(224.491,'J/(mol*K)'), comment="""Thermo group additivity estimation: group(Cs-(Cds-O2d)CsCsH) + group(Cs-CsHHH) + group(Cds-CdsCsH) + group(Cds-OdCsH) + group(Cds-CdsHH) + group(Cdd-CdsCds) + radical(C=C=CJ) + radical(C=CCJ(C)C=O) + radical(CJ3)"""),
)

species(
    label = 'C#C[CH][C]1[CH]C1[O](26155)',
    structure = SMILES('C#C[CH][C]1[CH]C1[O]'),
    E0 = (812.153,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
    molecularWeight = (92.0954,'amu'),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[1.85067,0.032439,2.56165e-05,-6.71564e-08,3.18006e-11,97770.7,24.0956], Tmin=(100,'K'), Tmax=(905.099,'K')), NASAPolynomial(coeffs=[16.1857,0.00849499,-1.17689e-08,-1.74081e-10,1.19401e-14,93561.6,-52.5503], Tmin=(905.099,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(812.153,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(249.434,'J/(mol*K)'), comment="""Thermo group additivity estimation: group(O2s-CsH) + group(Cs-CsCsCsH) + group(Cs-CsCsOsH) + group(Cs-CsCsHH) + group(Cs-CtCsHH) + group(Ct-CtCs) + group(Ct-CtH) + ring(Cyclopropane) + radical(CCJ(C)CO) + radical(CC(C)OJ) + radical(CCJCO) + radical(Sec_Propargyl)"""),
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
    label = '[CH]=C=C[C]([CH2])[C]=O(26156)',
    structure = SMILES('[CH]=C=C[C]([CH2])[C]=O'),
    E0 = (653.959,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
    molecularWeight = (92.0954,'amu'),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[1.33152,0.0567317,-5.83947e-05,3.19652e-08,-7.00391e-12,78750.8,26.2729], Tmin=(100,'K'), Tmax=(1107.59,'K')), NASAPolynomial(coeffs=[11.7365,0.0191549,-7.50504e-06,1.33443e-09,-9.01222e-14,76445.9,-24.9884], Tmin=(1107.59,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(653.959,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(245.277,'J/(mol*K)'), comment="""Thermo group additivity estimation: group(Cs-(Cds-O2d)CsCsH) + group(Cs-CsHHH) + group(Cds-CdsCsH) + group(Cds-OdCsH) + group(Cds-CdsHH) + group(Cdd-CdsCds) + radical(C=C=CJ) + radical(CJC(C)C=O) + radical(CC(C)CJ=O) + radical(C=CCJ(C)C=O)"""),
)

species(
    label = '[CH]=[C]C=C1[CH]O[CH]1(26157)',
    structure = SMILES('[CH]=[C]C=C1[CH]O[CH]1'),
    E0 = (729.302,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
    molecularWeight = (92.0954,'amu'),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[1.91217,0.0390122,-1.29762e-05,-6.64459e-09,4.086e-12,87795.7,22.9801], Tmin=(100,'K'), Tmax=(1115.04,'K')), NASAPolynomial(coeffs=[9.6875,0.0239351,-9.93387e-06,1.84403e-09,-1.28248e-13,85265.1,-18.9504], Tmin=(1115.04,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(729.302,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(253.591,'J/(mol*K)'), comment="""Thermo group additivity estimation: group(O2s-CsCs) + group(Cs-(Cds-Cds)OsHH) + group(Cs-(Cds-Cds)OsHH) + group(Cds-CdsCsCs) + group(Cds-Cds(Cds-Cds)H) + group(Cds-Cds(Cds-Cds)H) + group(Cds-CdsHH) + ring(Cyclobutane) + radical(C=CJC=C) + radical(C=CCJ(O)C) + radical(C=CCJ(O)C) + radical(Cds_P)"""),
)

species(
    label = '[CH]=C1[CH]C(C=O)=C1(26101)',
    structure = SMILES('[CH]=C1[CH]C(C=O)=C1'),
    E0 = (535.917,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
    molecularWeight = (92.0954,'amu'),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[1.67285,0.0363647,1.20412e-05,-4.6847e-08,2.13279e-11,64552.9,21.6982], Tmin=(100,'K'), Tmax=(985.418,'K')), NASAPolynomial(coeffs=[16.4741,0.0117857,-4.58615e-06,9.62677e-10,-7.68253e-14,59912.1,-58.2383], Tmin=(985.418,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(535.917,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(253.591,'J/(mol*K)'), comment="""Thermo group additivity estimation: group(Cs-(Cds-Cds)(Cds-Cds)HH) + group(Cd-CdCs(CO)) + group(Cds-Cds(Cds-Cds)Cs) + group(Cds-Cds(Cds-Cds)H) + group(Cds-O2d(Cds-Cds)H) + group(Cds-CdsHH) + ring(Cyclobutene) + radical(Cds_P) + radical(CCJCC=O)"""),
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
    collisionModel = TransportData(shapeIndex=2, epsilon=(3451.63,'J/mol'), sigma=(5.81624,'angstroms'), dipoleMoment=(0,'C*m'), polarizability=(0,'angstroms^3'), rotrelaxcollnum=0, comment="""Epsilon & sigma estimated with Tc=539.14 K, Pc=39.81 bar (from Joback method)"""),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[0.736887,0.0662785,-6.90866e-05,3.09145e-08,-3.59339e-12,98142.6,26.104], Tmin=(100,'K'), Tmax=(920.956,'K')), NASAPolynomial(coeffs=[16.833,0.011062,-3.08555e-06,4.61593e-10,-2.95335e-14,94554.7,-53.608], Tmin=(920.956,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(814.984,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(245.277,'J/(mol*K)'), comment="""Thermo group additivity estimation: group(O2s-Cs(Cds-Cd)) + group(Cs-(Cds-Cds)OsHH) + group(Cds-CdsCsH) + group(Cds-CdsOsH) + group(Cds-CdsHH) + group(Cds-CdsHH) + group(Cdd-CdsCds) + radical(Cds_S) + radical(C=C=CJ) + radical(C=CCJ(O)C) + radical(Cds_P)"""),
)

species(
    label = '[C]#C[CH]C([CH])C=O(26158)',
    structure = SMILES('[C]#C[CH]C([CH])C=O'),
    E0 = (998.722,'kJ/mol'),
    modes = [
        HarmonicOscillator(frequencies=([3025,407.5,1350,352.5,1380,1390,370,380,2900,435,2175,525,2782.5,750,1395,475,1775,1000,200,800,1066.67,1333.33,1600],'cm^-1')),
        HinderedRotor(inertia=(0.156089,'amu*angstrom^2'), symmetry=1, barrier=(3.5888,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.156089,'amu*angstrom^2'), symmetry=1, barrier=(3.5888,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.156089,'amu*angstrom^2'), symmetry=1, barrier=(3.5888,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.156089,'amu*angstrom^2'), symmetry=1, barrier=(3.5888,'kJ/mol'), semiclassical=False),
    ],
    spinMultiplicity = 5,
    opticalIsomers = 1,
    molecularWeight = (92.0954,'amu'),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[1.16546,0.0652686,-9.16784e-05,7.03185e-08,-2.10676e-11,120218,28.7206], Tmin=(100,'K'), Tmax=(927.49,'K')), NASAPolynomial(coeffs=[9.61448,0.0211984,-8.06212e-06,1.3444e-09,-8.46562e-14,118979,-9.63545], Tmin=(927.49,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(998.722,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(241.12,'J/(mol*K)'), comment="""Thermo group additivity estimation: group(Cs-(Cds-O2d)CsCsH) + group(Cs-CtCsHH) + group(Cs-CsHHH) + group(Cds-OdCsH) + group(Ct-CtCs) + group(Ct-CtH) + radical(CCJCC=O) + radical(CCJ2_triplet) + radical(Acetyl)"""),
)

species(
    label = '[C]#C[CH]C([CH2])=C[O](26159)',
    structure = SMILES('[C]#C[CH]C([CH2])=C[O]'),
    E0 = (834.93,'kJ/mol'),
    modes = [
        HarmonicOscillator(frequencies=([3025,407.5,1350,352.5,3010,987.5,1337.5,450,1655,350,440,435,1725,2175,525,3000,3100,440,815,1455,1000,412.389,413.465,413.556],'cm^-1')),
        HinderedRotor(inertia=(0.000995216,'amu*angstrom^2'), symmetry=1, barrier=(0.119627,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.37081,'amu*angstrom^2'), symmetry=1, barrier=(44.4608,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.367074,'amu*angstrom^2'), symmetry=1, barrier=(44.4992,'kJ/mol'), semiclassical=False),
    ],
    spinMultiplicity = 5,
    opticalIsomers = 1,
    molecularWeight = (92.0954,'amu'),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[0.807054,0.0621114,-6.58844e-05,3.49674e-08,-7.14751e-12,100541,27.0555], Tmin=(100,'K'), Tmax=(1261.97,'K')), NASAPolynomial(coeffs=[15.9393,0.0118326,-3.37056e-06,4.89389e-10,-2.93501e-14,96905.8,-48.7392], Tmin=(1261.97,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(834.93,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(245.277,'J/(mol*K)'), comment="""Thermo group additivity estimation: group(O2s-(Cds-Cd)H) + group(Cs-(Cds-Cds)CtHH) + group(Cs-(Cds-Cds)HHH) + group(Cds-CdsCsCs) + group(Cds-CdsOsH) + group(Ct-CtCs) + group(Ct-CtH) + radical(Acetyl) + radical(C=COJ) + radical(Allyl_P) + radical(Cs_S)"""),
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
    label = '[CH]=C1C([O])C1C#C(26161)',
    structure = SMILES('[CH]=C1C([O])C1C#C'),
    E0 = (704.847,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
    molecularWeight = (92.0954,'amu'),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[0.98612,0.0571845,-4.34801e-05,5.90522e-09,4.42314e-12,84890.2,20.7503], Tmin=(100,'K'), Tmax=(962.681,'K')), NASAPolynomial(coeffs=[17.2464,0.0104965,-3.25853e-06,5.75435e-10,-4.20788e-14,80792.3,-62.1013], Tmin=(962.681,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(704.847,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(253.591,'J/(mol*K)'), comment="""Thermo group additivity estimation: group(O2s-CsH) + group(Cs-(Cds-Cds)CtCsH) + group(Cs-(Cds-Cds)CsOsH) + group(Cds-CdsCsCs) + group(Cds-CdsHH) + group(Ct-CtCs) + group(Ct-CtH) + ring(Methylene_cyclopropane) + radical(Cds_P) + radical(CC(C)OJ)"""),
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
    label = '[CH]=[C]C([C]=[CH])C=O(26163)',
    structure = SMILES('[CH]=[C]C([C]=[CH])C=O'),
    E0 = (942.746,'kJ/mol'),
    modes = [
        HarmonicOscillator(frequencies=([3115,3125,620,680,785,800,1600,1700,2782.5,750,1395,475,1775,1000,1380,1390,370,380,2900,435,1670,1700,300,440],'cm^-1')),
        HinderedRotor(inertia=(0.460331,'amu*angstrom^2'), symmetry=1, barrier=(10.5839,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.459878,'amu*angstrom^2'), symmetry=1, barrier=(10.5735,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.460207,'amu*angstrom^2'), symmetry=1, barrier=(10.5811,'kJ/mol'), semiclassical=False),
    ],
    spinMultiplicity = 5,
    opticalIsomers = 1,
    molecularWeight = (92.0954,'amu'),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[1.0278,0.0758735,-0.000135167,1.27075e-07,-4.52654e-11,113483,27.3265], Tmin=(100,'K'), Tmax=(863.054,'K')), NASAPolynomial(coeffs=[5.87708,0.0300416,-1.49162e-05,2.83001e-09,-1.92135e-13,113516,9.68535], Tmin=(863.054,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(942.746,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(245.277,'J/(mol*K)'), comment="""Thermo group additivity estimation: group(Cs-(Cds-O2d)CsCsH) + group(Cds-CdsCsH) + group(Cds-CdsCsH) + group(Cds-OdCsH) + group(Cds-CdsHH) + group(Cds-CdsHH) + radical(Cds_S) + radical(Cds_P) + radical(Cds_P) + radical(Cds_S)"""),
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
    E0 = (654.74,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS2',
    E0 = (879.316,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS3',
    E0 = (1134.69,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS4',
    E0 = (852.89,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS5',
    E0 = (915.801,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS6',
    E0 = (796.136,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS7',
    E0 = (841.152,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS8',
    E0 = (772.179,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS9',
    E0 = (1039,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS10',
    E0 = (1104.39,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS11',
    E0 = (1103.47,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS12',
    E0 = (878.408,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS13',
    E0 = (843.375,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS14',
    E0 = (769.607,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS15',
    E0 = (729.257,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS16',
    E0 = (679.714,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS17',
    E0 = (1201.88,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS18',
    E0 = (1141.13,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS19',
    E0 = (956.24,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS20',
    E0 = (1343.88,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS21',
    E0 = (1265.92,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS22',
    E0 = (945.783,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS23',
    E0 = (974.729,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS24',
    E0 = (1143.75,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS25',
    E0 = (942.058,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS26',
    E0 = (905.13,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS27',
    E0 = (685.773,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS28',
    E0 = (743.709,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS29',
    E0 = (662.272,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS30',
    E0 = (1227.36,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS31',
    E0 = (899.116,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS32',
    E0 = (1028.8,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS33',
    E0 = (930.161,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS34',
    E0 = (910.942,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS35',
    E0 = (732.988,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS36',
    E0 = (1180.19,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS37',
    E0 = (850.76,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS38',
    E0 = (881.92,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS39',
    E0 = (829.652,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS40',
    E0 = (1072.27,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS41',
    E0 = (1170.58,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS42',
    E0 = (1014.69,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS43',
    E0 = (1473.06,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS44',
    E0 = (1374.76,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS45',
    E0 = (1022.59,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS46',
    E0 = (1024.47,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS47',
    E0 = (1340.49,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS48',
    E0 = (1295.7,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS49',
    E0 = (654.74,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS50',
    E0 = (1294.23,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS51',
    E0 = (817.367,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS52',
    E0 = (1013.5,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS53',
    E0 = (797.58,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS54',
    E0 = (785.654,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS55',
    E0 = (792.796,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS56',
    E0 = (1043.1,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS57',
    E0 = (657.669,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS58',
    E0 = (1166.04,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS59',
    E0 = (812.153,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS60',
    E0 = (975.117,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS61',
    E0 = (897.507,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS62',
    E0 = (846.869,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS63',
    E0 = (1131.62,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS64',
    E0 = (792.796,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS65',
    E0 = (663.025,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS66',
    E0 = (1128.78,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS67',
    E0 = (1151.2,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS68',
    E0 = (944.147,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS69',
    E0 = (705.367,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS70',
    E0 = (704.847,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS71',
    E0 = (662.272,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS72',
    E0 = (1065.5,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

reaction(
    label = 'reaction1',
    reactants = ['[CH]C(C=C=[CH])=C[O](22639)'],
    products = ['C#CC=O(21959)', '[CH]=C=[CH](18734)'],
    transitionState = 'TS1',
    kinetics = Arrhenius(A=(5e+12,'s^-1'), n=0, Ea=(0,'kJ/mol'), T0=(1,'K'), Tmin=(300,'K'), Tmax=(1500,'K'), comment="""Exact match found for rate rule [RJJ]
Euclidian distance = 0
family: 1,4_Linear_birad_scission"""),
)

reaction(
    label = 'reaction2',
    reactants = ['[CH]C(C=C=[CH])=C[O](22639)'],
    products = ['[CH]C1=COC1[C]=[CH](25320)'],
    transitionState = 'TS2',
    kinetics = Arrhenius(A=(1e+10,'s^-1'), n=0, Ea=(224.575,'kJ/mol'), T0=(1,'K'), comment="""Estimated using template [R5_SD_D;doublebond_intra;radadd_intra] for rate rule [R5_SD_D;doublebond_intra;radadd_intra_O]
Euclidian distance = 1.0
family: Intra_R_Add_Exocyclic
Ea raised from 222.6 to 224.6 kJ/mol to match endothermicity of reaction."""),
)

reaction(
    label = 'reaction3',
    reactants = ['[CH]C(C=C=[CH])=C[O](22639)'],
    products = ['[CH]C1([CH][O])C=C=C1(26121)'],
    transitionState = 'TS3',
    kinetics = Arrhenius(A=(5.41e+10,'s^-1'), n=0.21, Ea=(479.947,'kJ/mol'), T0=(1,'K'), comment="""Estimated using template [R5;doublebond_intra;radadd_intra_cdsingleH] for rate rule [R5_MM;doublebond_intra;radadd_intra_cdsingleH]
Euclidian distance = 1.0
family: Intra_R_Add_Exocyclic
Ea raised from 479.3 to 479.9 kJ/mol to match endothermicity of reaction."""),
)

reaction(
    label = 'reaction4',
    reactants = ['H(8)', '[CH]C([C]=O)=CC#C(26122)'],
    products = ['[CH]C(C=C=[CH])=C[O](22639)'],
    transitionState = 'TS4',
    kinetics = Arrhenius(A=(7048.5,'m^3/(mol*s)'), n=1.43107, Ea=(2.92155,'kJ/mol'), T0=(1,'K'), Tmin=(303.03,'K'), Tmax=(2000,'K'), comment="""From training reaction 188 used for Ck_O;HJ
Exact match found for rate rule [Ck_O;HJ]
Euclidian distance = 0
family: R_Addition_MultipleBond"""),
)

reaction(
    label = 'reaction5',
    reactants = ['[CH]C(=[C]O)C=C=[CH](26123)'],
    products = ['[CH]C(C=C=[CH])=C[O](22639)'],
    transitionState = 'TS5',
    kinetics = Arrhenius(A=(4.96519e+09,'s^-1'), n=1.05826, Ea=(162.779,'kJ/mol'), T0=(1,'K'), comment="""Estimated using template [R2H_S;Cd_rad_out_Cd;XH_out] for rate rule [R2H_S;Cd_rad_out_Cd;O_H_out]
Euclidian distance = 1.0
family: intra_H_migration"""),
)

reaction(
    label = 'reaction6',
    reactants = ['[CH]C(C=C=[CH])=C[O](22639)'],
    products = ['[CH]C(C#C[CH2])=C[O](24938)'],
    transitionState = 'TS6',
    kinetics = Arrhenius(A=(1.59701e+06,'s^-1'), n=1.91493, Ea=(141.395,'kJ/mol'), T0=(1,'K'), comment="""Estimated using average of templates [RnH;Cd_rad_out_singleH;Cd_H_out_singleDe] + [R3H;Cd_rad_out_singleH;XH_out] for rate rule [R3H;Cd_rad_out_singleH;Cd_H_out_singleDe]
Euclidian distance = 2.0
family: intra_H_migration"""),
)

reaction(
    label = 'reaction7',
    reactants = ['[CH]C(=C=C=[CH])[CH]O(26124)'],
    products = ['[CH]C(C=C=[CH])=C[O](22639)'],
    transitionState = 'TS7',
    kinetics = Arrhenius(A=(1.11e+08,'s^-1'), n=1.1915, Ea=(103.605,'kJ/mol'), T0=(1,'K'), comment="""Estimated using template [R4H_SDS;Cd_rad_out_double;XH_out] for rate rule [R4H_SDS;Cd_rad_out_double;O_H_out]
Euclidian distance = 1.0
family: intra_H_migration"""),
)

reaction(
    label = 'reaction8',
    reactants = ['[CH]C([C]=O)=C[C]=C(24939)'],
    products = ['[CH]C(C=C=[CH])=C[O](22639)'],
    transitionState = 'TS8',
    kinetics = Arrhenius(A=(33613.8,'s^-1'), n=2.10442, Ea=(111.806,'kJ/mol'), T0=(1,'K'), comment="""Estimated using template [R5H;Y_rad_out;Cd_H_out_singleH] for rate rule [R5H;Cd_rad_out;Cd_H_out_singleH]
Euclidian distance = 1.0
Multiplied by reaction path degeneracy 2.0
family: intra_H_migration"""),
)

reaction(
    label = 'reaction9',
    reactants = ['[CH][C]=C[O](21209)', '[CH]=C=[CH](18734)'],
    products = ['[CH]C(C=C=[CH])=C[O](22639)'],
    transitionState = 'TS9',
    kinetics = Arrhenius(A=(81155.9,'m^3/(mol*s)'), n=0.702818, Ea=(0,'kJ/mol'), T0=(1,'K'), comment="""Estimated using an average for rate rule [Cd_allenic;Y_rad]
Euclidian distance = 0
Multiplied by reaction path degeneracy 2.0
family: R_Recombination
Ea raised from -11.4 to 0 kJ/mol."""),
)

reaction(
    label = 'reaction10',
    reactants = ['H(8)', '[CH]C([C]=C=[CH])=C[O](26125)'],
    products = ['[CH]C(C=C=[CH])=C[O](22639)'],
    transitionState = 'TS10',
    kinetics = Arrhenius(A=(2.048e+13,'cm^3/(mol*s)'), n=0.206, Ea=(0,'kJ/mol'), T0=(1,'K'), comment="""From training reaction 44 used for Cd_allenic;H_rad
Exact match found for rate rule [Cd_allenic;H_rad]
Euclidian distance = 0
family: R_Recombination
Ea raised from -0.7 to 0 kJ/mol."""),
)

reaction(
    label = 'reaction11',
    reactants = ['H(8)', '[CH][C]([C]=O)C=C=[CH](26126)'],
    products = ['[CH]C(C=C=[CH])=C[O](22639)'],
    transitionState = 'TS11',
    kinetics = Arrhenius(A=(4.34078e+06,'m^3/(mol*s)'), n=0.278577, Ea=(0,'kJ/mol'), T0=(1,'K'), comment="""Estimated using an average for rate rule [Y_rad;H_rad]
Euclidian distance = 0
family: R_Recombination
Ea raised from -1.4 to 0 kJ/mol."""),
)

reaction(
    label = 'reaction12',
    reactants = ['[CH]C(C=C=[CH])=C[O](22639)'],
    products = ['[CH]C1([CH]O1)C=C=[CH](26127)'],
    transitionState = 'TS12',
    kinetics = Arrhenius(A=(2.42595e+10,'s^-1'), n=0.7335, Ea=(223.668,'kJ/mol'), T0=(1,'K'), comment="""Estimated using template [R3_D;doublebond_intra_pri;radadd_intra] for rate rule [R3_D;doublebond_intra_pri;radadd_intra_O]
Euclidian distance = 1.0
family: Intra_R_Add_Endocyclic
Ea raised from 223.0 to 223.7 kJ/mol to match endothermicity of reaction."""),
)

reaction(
    label = 'reaction13',
    reactants = ['[CH]C(C=C=[CH])=C[O](22639)'],
    products = ['[CH]C(=C[O])C1[C]=C1(26128)'],
    transitionState = 'TS13',
    kinetics = Arrhenius(A=(4.2354e+12,'s^-1'), n=-0.1205, Ea=(188.635,'kJ/mol'), T0=(1,'K'), comment="""Estimated using average of templates [Rn;doublebond_intra_CdCdd;radadd_intra_cdsingleH] + [R3;doublebond_intra_CdCdd;radadd_intra_cdsingle] for rate rule [R3;doublebond_intra_CdCdd;radadd_intra_cdsingleH]
Euclidian distance = 1.0
family: Intra_R_Add_Endocyclic
Ea raised from 188.2 to 188.6 kJ/mol to match endothermicity of reaction."""),
)

reaction(
    label = 'reaction14',
    reactants = ['[CH]C(C=C=[CH])=C[O](22639)'],
    products = ['[CH]C1=COC([CH])=C1(26106)'],
    transitionState = 'TS14',
    kinetics = Arrhenius(A=(8.9814e+11,'s^-1'), n=0.0209575, Ea=(114.866,'kJ/mol'), T0=(1,'K'), comment="""Estimated using template [R5_SD_D;doublebond_intra;radadd_intra] for rate rule [R5_SD_D;doublebond_intra;radadd_intra_O]
Euclidian distance = 1.0
family: Intra_R_Add_Endocyclic"""),
)

reaction(
    label = 'reaction15',
    reactants = ['[CH]C(C=C=[CH])=C[O](22639)'],
    products = ['[CH]=C1C=[C][CH]C1[O](26129)'],
    transitionState = 'TS15',
    kinetics = Arrhenius(A=(9.13977e+10,'s^-1'), n=0.0396934, Ea=(74.5166,'kJ/mol'), T0=(1,'K'), comment="""Estimated using template [R5_linear;doublebond_intra;radadd_intra_cdsingleH] for rate rule [R5_MM;doublebond_intra;radadd_intra_cdsingleH]
Euclidian distance = 1.0
family: Intra_R_Add_Endocyclic
Ea raised from 68.6 to 74.5 kJ/mol to match endothermicity of reaction."""),
)

reaction(
    label = 'reaction16',
    reactants = ['[CH]C(C=C=[CH])=C[O](22639)'],
    products = ['[CH]C(=C=O)C=C=C(24945)'],
    transitionState = 'TS16',
    kinetics = Arrhenius(A=(2.1261e+09,'s^-1'), n=0.137, Ea=(24.9733,'kJ/mol'), T0=(1,'K'), comment="""Estimated using template [R5;Y_rad;XH_Rrad] for rate rule [R5radExo;Y_rad;XH_Rrad]
Euclidian distance = 1.0
family: Intra_Disproportionation"""),
)

reaction(
    label = 'reaction17',
    reactants = ['[CH]C([C]=C=[CH])[CH][O](26130)'],
    products = ['[CH]C(C=C=[CH])=C[O](22639)'],
    transitionState = 'TS17',
    kinetics = Arrhenius(A=(1.949e+11,'s^-1'), n=0.486, Ea=(22.8614,'kJ/mol'), T0=(1,'K'), comment="""Estimated using an average for rate rule [R2radExo;Y_rad_De;XH_Rrad]
Euclidian distance = 0
family: Intra_Disproportionation"""),
)

reaction(
    label = 'reaction18',
    reactants = ['[CH]C(=[C][O])C[C]=[CH](26131)'],
    products = ['[CH]C(C=C=[CH])=C[O](22639)'],
    transitionState = 'TS18',
    kinetics = Arrhenius(A=(1.4874e+09,'s^-1'), n=1.045, Ea=(63.4002,'kJ/mol'), T0=(1,'K'), comment="""Estimated using an average for rate rule [R3radExo;Y_rad;XH_Rrad]
Euclidian distance = 0
Multiplied by reaction path degeneracy 2.0
family: Intra_Disproportionation"""),
)

reaction(
    label = 'reaction19',
    reactants = ['[CH]C=C[C]([CH])[C]=O(26132)'],
    products = ['[CH]C(C=C=[CH])=C[O](22639)'],
    transitionState = 'TS19',
    kinetics = Arrhenius(A=(5.14222e+08,'s^-1'), n=0.311, Ea=(24.9733,'kJ/mol'), T0=(1,'K'), comment="""Estimated using template [R4;Y_rad;XH_Rrad] for rate rule [R4radEndo;Y_rad;XH_Rrad]
Euclidian distance = 1.0
family: Intra_Disproportionation"""),
)

reaction(
    label = 'reaction20',
    reactants = ['O(T)(63)', '[CH]C(=[CH])C=C=[CH](21251)'],
    products = ['[CH]C(C=C=[CH])=C[O](22639)'],
    transitionState = 'TS20',
    kinetics = Arrhenius(A=(187219,'m^3/(mol*s)'), n=1.13083, Ea=(163.847,'kJ/mol'), T0=(1,'K'), Tmin=(303.03,'K'), Tmax=(2000,'K'), comment="""Estimated using template [Y_rad;O_birad] for rate rule [Cd_pri_rad;O_birad]
Euclidian distance = 2.0
Multiplied by reaction path degeneracy 2.0
family: Birad_R_Recombination"""),
)

reaction(
    label = 'reaction21',
    reactants = ['H(8)', '[C]#C[CH]C([CH])=C[O](26133)'],
    products = ['[CH]C(C=C=[CH])=C[O](22639)'],
    transitionState = 'TS21',
    kinetics = Arrhenius(A=(1e+07,'m^3/(mol*s)'), n=0, Ea=(0,'kJ/mol'), T0=(1,'K'), comment="""Estimated using an average for rate rule [H_rad;Birad]
Euclidian distance = 0
family: Birad_R_Recombination"""),
)

reaction(
    label = 'reaction22',
    reactants = ['[CH](2815)', 'C#CC=C=C[O](23376)'],
    products = ['[CH]C(C=C=[CH])=C[O](22639)'],
    transitionState = 'TS22',
    kinetics = Arrhenius(A=(1149,'m^3/(mol*s)'), n=1.595, Ea=(16.7946,'kJ/mol'), T0=(1,'K'), comment="""Estimated using template [Ca_Cds-OneDeH;YJ] for rate rule [Ca_Cds-CtH;CH_quartet]
Euclidian distance = 2.2360679775
family: R_Addition_MultipleBond"""),
)

reaction(
    label = 'reaction23',
    reactants = ['[CH]C(=C=O)C[C]=[CH](24774)'],
    products = ['[CH]C(C=C=[CH])=C[O](22639)'],
    transitionState = 'TS23',
    kinetics = Arrhenius(A=(1.47715e+10,'s^-1'), n=0.8, Ea=(147.277,'kJ/mol'), T0=(1,'K'), comment="""Estimated using template [R3H_DS;Cd_rad_out;Cs_H_out_1H] for rate rule [R3H_DS;Cd_rad_out;Cs_H_out_H/Ct]
Euclidian distance = 2.0
Multiplied by reaction path degeneracy 2.0
family: intra_H_migration"""),
)

reaction(
    label = 'reaction24',
    reactants = ['[C]=[C]CC(=[CH])C=O(24783)'],
    products = ['[CH]C(C=C=[CH])=C[O](22639)'],
    transitionState = 'TS24',
    kinetics = Arrhenius(A=(746634,'s^-1'), n=2.07188, Ea=(142.765,'kJ/mol'), T0=(1,'K'), comment="""Estimated using template [R3H;Y_rad_out;Cs_H_out_H/Cd] for rate rule [R3H_TS;Ct_rad_out;Cs_H_out_H/Cd]
Euclidian distance = 2.2360679775
Multiplied by reaction path degeneracy 2.0
family: intra_H_migration"""),
)

reaction(
    label = 'reaction25',
    reactants = ['[C]#CC=C([CH])[CH]O(26134)'],
    products = ['[CH]C(C=C=[CH])=C[O](22639)'],
    transitionState = 'TS25',
    kinetics = Arrhenius(A=(19101.7,'s^-1'), n=1.79759, Ea=(74.4436,'kJ/mol'), T0=(1,'K'), comment="""Estimated using template [R6Hall;Y_rad_out;XH_out] for rate rule [R6HJ_2;Ct_rad_out;O_H_out]
Euclidian distance = 1.73205080757
family: intra_H_migration"""),
)

reaction(
    label = 'reaction26',
    reactants = ['[CH]C(C=C=[CH])=C[O](22639)'],
    products = ['[CH][C]1C([O])C1C#C(26135)'],
    transitionState = 'TS26',
    kinetics = Arrhenius(A=(3.473e+12,'s^-1'), n=0.247, Ea=(250.39,'kJ/mol'), T0=(1,'K'), comment="""Estimated using template [R3_D;doublebond_intra;radadd_intra_csHDe] for rate rule [R3_D;doublebond_intra;radadd_intra_csHCt]
Euclidian distance = 1.0
family: Intra_R_Add_Endocyclic
Ea raised from 246.4 to 250.4 kJ/mol to match endothermicity of reaction."""),
)

reaction(
    label = 'reaction27',
    reactants = ['[CH]C(C=C=[CH])=C[O](22639)'],
    products = ['[CH]C1[CH]OC=[C]C=1(25326)'],
    transitionState = 'TS27',
    kinetics = Arrhenius(A=(9.23539e+09,'s^-1'), n=0.445806, Ea=(31.0324,'kJ/mol'), T0=(1,'K'), comment="""Estimated using template [R6_linear;triplebond_intra_H;radadd_intra] for rate rule [R6_linear;triplebond_intra_H;radadd_intra_O]
Euclidian distance = 1.0
family: Intra_R_Add_Endocyclic"""),
)

reaction(
    label = 'reaction28',
    reactants = ['[CH]C(C=C=[CH])=C[O](22639)'],
    products = ['[CH]C(=C=O)CC#C(26136)'],
    transitionState = 'TS28',
    kinetics = Arrhenius(A=(2.6374e+09,'s^-1'), n=0.37, Ea=(88.9686,'kJ/mol'), T0=(1,'K'), comment="""Estimated using average of templates [R3;Y_rad_De;XH_Rrad] + [R3radExo;Y_rad;XH_Rrad] for rate rule [R3radExo;Y_rad_De;XH_Rrad]
Euclidian distance = 1.0
family: Intra_Disproportionation"""),
)

reaction(
    label = 'reaction29',
    reactants = ['[CH]C(C=C=[CH])=C[O](22639)'],
    products = ['[CH]C1=COC1C#C(25333)'],
    transitionState = 'TS29',
    kinetics = Arrhenius(A=(2e+12,'s^-1'), n=0, Ea=(7.5312,'kJ/mol'), T0=(1,'K'), comment="""Estimated using template [Rn;C_rad_out_H/OneDe;Ypri_rad_out] for rate rule [R4_SDS;C_rad_out_H/OneDe;Opri_rad]
Euclidian distance = 2.2360679775
family: Birad_recombination"""),
)

reaction(
    label = 'reaction30',
    reactants = ['[C]#C(5143)', '[CH]C([CH])=C[O](21383)'],
    products = ['[CH]C(C=C=[CH])=C[O](22639)'],
    transitionState = 'TS30',
    kinetics = Arrhenius(A=(2.29708e+06,'m^3/(mol*s)'), n=0.575199, Ea=(34.3157,'kJ/mol'), T0=(1,'K'), comment="""Estimated using template [Y_rad;Birad] for rate rule [Ct_rad;Birad]
Euclidian distance = 1.0
Multiplied by reaction path degeneracy 2.0
family: Birad_R_Recombination"""),
)

reaction(
    label = 'reaction31',
    reactants = ['H(8)', '[CH]C(=C=C=[CH])C=O(26137)'],
    products = ['[CH]C(C=C=[CH])=C[O](22639)'],
    transitionState = 'TS31',
    kinetics = Arrhenius(A=(8.22e+08,'cm^3/(mol*s)'), n=1.533, Ea=(7.77387,'kJ/mol'), T0=(1,'K'), comment="""From training reaction 192 used for Cd_R;HJ
Exact match found for rate rule [Cd_R;HJ]
Euclidian distance = 0
family: R_Addition_MultipleBond"""),
)

reaction(
    label = 'reaction32',
    reactants = ['[CH]C([C]=C=[CH])C=O(26138)'],
    products = ['[CH]C(C=C=[CH])=C[O](22639)'],
    transitionState = 'TS32',
    kinetics = Arrhenius(A=(3.9024e+09,'s^-1'), n=1.09833, Ea=(177.686,'kJ/mol'), T0=(1,'K'), comment="""Estimated using an average for rate rule [R2H_S;Cd_rad_out_double;XH_out]
Euclidian distance = 0
family: intra_H_migration"""),
)

reaction(
    label = 'reaction33',
    reactants = ['[CH]C([C]=O)C=C=[CH](26139)'],
    products = ['[CH]C(C=C=[CH])=C[O](22639)'],
    transitionState = 'TS33',
    kinetics = Arrhenius(A=(4.89315e+07,'s^-1'), n=1.5961, Ea=(158.199,'kJ/mol'), T0=(1,'K'), comment="""Estimated using template [R2H_S;Y_rad_out;XH_out] for rate rule [R2H_S;CO_rad_out;XH_out]
Euclidian distance = 1.0
family: intra_H_migration"""),
)

reaction(
    label = 'reaction34',
    reactants = ['[CH]C(C=C=[CH])=C[O](22639)'],
    products = ['[CH]C1=CC1([CH])C=O(26140)'],
    transitionState = 'TS34',
    kinetics = Arrhenius(A=(4.00692e+11,'s^-1'), n=0.347401, Ea=(256.201,'kJ/mol'), T0=(1,'K'), comment="""Estimated using an average for rate rule [R3_D;doublebond_intra;radadd_intra_cs]
Euclidian distance = 0
family: Intra_R_Add_Endocyclic
Ea raised from 254.2 to 256.2 kJ/mol to match endothermicity of reaction."""),
)

reaction(
    label = 'reaction35',
    reactants = ['[CH]C(C=C=[CH])=C[O](22639)'],
    products = ['[CH]=C(C#C[CH2])C=O(24935)'],
    transitionState = 'TS35',
    kinetics = Arrhenius(A=(2.00399e+09,'s^-1'), n=0.37, Ea=(78.2471,'kJ/mol'), T0=(1,'K'), comment="""Estimated using average of templates [R3;Y_rad;XH_Rrad_De] + [R3radExo;Y_rad;XH_Rrad] for rate rule [R3radExo;Y_rad;XH_Rrad_De]
Euclidian distance = 1.0
family: Intra_Disproportionation"""),
)

reaction(
    label = 'reaction36',
    reactants = ['[CH]=O(373)', '[CH][C]=C[C]=[CH](21426)'],
    products = ['[CH]C(C=C=[CH])=C[O](22639)'],
    transitionState = 'TS36',
    kinetics = Arrhenius(A=(2.29708e+06,'m^3/(mol*s)'), n=0.575199, Ea=(34.3157,'kJ/mol'), T0=(1,'K'), comment="""Estimated using template [Y_rad;Birad] for rate rule [CO_pri_rad;Birad]
Euclidian distance = 2.0
Multiplied by reaction path degeneracy 2.0
family: Birad_R_Recombination"""),
)

reaction(
    label = 'reaction37',
    reactants = ['[CH]C(C=C=[CH])=C[O](22639)'],
    products = ['[CH]=C1[CH]C(=[CH])C1[O](26113)'],
    transitionState = 'TS37',
    kinetics = Arrhenius(A=(1e+10,'s^-1'), n=0, Ea=(196.02,'kJ/mol'), T0=(1,'K'), comment="""Estimated using template [R5_SD;multiplebond_intra;radadd_intra] for rate rule [R5_SD_CO;carbonylbond_intra_H;radadd_intra]
Euclidian distance = 2.2360679775
family: Intra_R_Add_Exocyclic"""),
)

reaction(
    label = 'reaction38',
    reactants = ['[CH]C=C=C([CH])C=O(26141)'],
    products = ['[CH]C(C=C=[CH])=C[O](22639)'],
    transitionState = 'TS38',
    kinetics = Arrhenius(A=(4.96519e+09,'s^-1'), n=1.05826, Ea=(162.779,'kJ/mol'), T0=(1,'K'), comment="""Estimated using an average for rate rule [R2H_S;Cd_rad_out_Cd;XH_out]
Euclidian distance = 0
family: intra_H_migration"""),
)

reaction(
    label = 'reaction39',
    reactants = ['[CH]C(=C=O)[CH]C=[CH](26142)'],
    products = ['[CH]C(C=C=[CH])=C[O](22639)'],
    transitionState = 'TS39',
    kinetics = Arrhenius(A=(8.08094e+07,'s^-1'), n=1.1965, Ea=(138.911,'kJ/mol'), T0=(1,'K'), comment="""Estimated using template [R4H_SDS;Y_rad_out;XH_out] for rate rule [R4H_SDS;CO_rad_out;XH_out]
Euclidian distance = 1.0
family: intra_H_migration"""),
)

reaction(
    label = 'reaction40',
    reactants = ['[CH][C]=CC([CH])[C]=O(26143)'],
    products = ['[CH]C(C=C=[CH])=C[O](22639)'],
    transitionState = 'TS40',
    kinetics = Arrhenius(A=(1.949e+11,'s^-1'), n=0.486, Ea=(22.8614,'kJ/mol'), T0=(1,'K'), comment="""Estimated using an average for rate rule [R2radExo;Y_rad;XH_Rrad]
Euclidian distance = 0
family: Intra_Disproportionation"""),
)

reaction(
    label = 'reaction41',
    reactants = ['[CH]C(=[C][C]=[CH])C[O](26144)'],
    products = ['[CH]C(C=C=[CH])=C[O](22639)'],
    transitionState = 'TS41',
    kinetics = Arrhenius(A=(1.4874e+09,'s^-1'), n=1.045, Ea=(63.4002,'kJ/mol'), T0=(1,'K'), comment="""Estimated using an average for rate rule [R3radExo;Y_rad;XH_Rrad]
Euclidian distance = 0
Multiplied by reaction path degeneracy 2.0
family: Intra_Disproportionation"""),
)

reaction(
    label = 'reaction42',
    reactants = ['[CH][C]=[C]C([CH])=CO(26145)'],
    products = ['[CH]C(C=C=[CH])=C[O](22639)'],
    transitionState = 'TS42',
    kinetics = Arrhenius(A=(5.14222e+08,'s^-1'), n=0.311, Ea=(24.9733,'kJ/mol'), T0=(1,'K'), comment="""Estimated using template [R4;Y_rad;XH_Rrad] for rate rule [R4radEndo;Y_rad;XH_Rrad]
Euclidian distance = 1.0
family: Intra_Disproportionation"""),
)

reaction(
    label = 'reaction43',
    reactants = ['[C-]#[O+](374)', '[CH][CH]C=C=[CH](19252)'],
    products = ['[CH]C(C=C=[CH])=C[O](22639)'],
    transitionState = 'TS43',
    kinetics = Arrhenius(A=(0.0591985,'m^3/(mol*s)'), n=2.3675, Ea=(305.306,'kJ/mol'), T0=(1,'K'), comment="""Estimated using an average for rate rule [CO;R_H]
Euclidian distance = 0
family: 1,2_Insertion_CO"""),
)

reaction(
    label = 'reaction44',
    reactants = ['[C]=[CH](18830)', '[CH]C(=[CH])C=O(21223)'],
    products = ['[CH]C(C=C=[CH])=C[O](22639)'],
    transitionState = 'TS44',
    kinetics = Arrhenius(A=(2.29708e+06,'m^3/(mol*s)'), n=0.575199, Ea=(34.3157,'kJ/mol'), T0=(1,'K'), comment="""Estimated using template [Y_rad;Birad] for rate rule [Cd_pri_rad;Birad]
Euclidian distance = 2.0
Multiplied by reaction path degeneracy 2.0
family: Birad_R_Recombination"""),
)

reaction(
    label = 'reaction45',
    reactants = ['[CH]C(=C=C=[CH])C[O](26146)'],
    products = ['[CH]C(C=C=[CH])=C[O](22639)'],
    transitionState = 'TS45',
    kinetics = Arrhenius(A=(1.65652e+07,'s^-1'), n=1.67955, Ea=(176.633,'kJ/mol'), T0=(1,'K'), comment="""Estimated using template [R3H_DS;Cd_rad_out_singleDe;XH_out] for rate rule [R3H_DS;Cd_rad_out_singleDe_Ct;XH_out]
Euclidian distance = 1.0
Multiplied by reaction path degeneracy 2.0
family: intra_H_migration"""),
)

reaction(
    label = 'reaction46',
    reactants = ['[C]#CC=C([CH])C[O](26147)'],
    products = ['[CH]C(C=C=[CH])=C[O](22639)'],
    transitionState = 'TS46',
    kinetics = Arrhenius(A=(1.2946e+07,'s^-1'), n=1.39981, Ea=(48.4448,'kJ/mol'), T0=(1,'K'), comment="""Estimated using template [R5H_RSMS;Y_rad_out;XH_out] for rate rule [R5H_TSMS;Ct_rad_out;XH_out]
Euclidian distance = 1.41421356237
Multiplied by reaction path degeneracy 2.0
family: intra_H_migration"""),
)

reaction(
    label = 'reaction47',
    reactants = ['[C]#CC[C]([CH])[CH][O](26148)'],
    products = ['[CH]C(C=C=[CH])=C[O](22639)'],
    transitionState = 'TS47',
    kinetics = Arrhenius(A=(1.4874e+09,'s^-1'), n=1.045, Ea=(63.4002,'kJ/mol'), T0=(1,'K'), comment="""Estimated using an average for rate rule [R3radExo;Y_rad;XH_Rrad]
Euclidian distance = 0
Multiplied by reaction path degeneracy 2.0
family: Intra_Disproportionation"""),
)

reaction(
    label = 'reaction48',
    reactants = ['[C]#C[CH]C([CH])[CH][O](26149)'],
    products = ['[CH]C(C=C=[CH])=C[O](22639)'],
    transitionState = 'TS48',
    kinetics = Arrhenius(A=(5.14222e+08,'s^-1'), n=0.311, Ea=(24.9733,'kJ/mol'), T0=(1,'K'), comment="""Estimated using template [R4;Y_rad;XH_Rrad] for rate rule [R4radEndo;Y_rad;XH_Rrad]
Euclidian distance = 1.0
family: Intra_Disproportionation"""),
)

reaction(
    label = 'reaction49',
    reactants = ['[CH]C(C=C=[CH])=C[O](22639)'],
    products = ['[CH]C(C=O)=CC#C(22634)'],
    transitionState = 'TS49',
    kinetics = Arrhenius(A=(1e+10,'s^-1'), n=0, Ea=(0,'kJ/mol'), T0=(1,'K'), comment="""Estimated using an average for rate rule [NOS]
Euclidian distance = 0
family: 1,2-Birad_to_alkene"""),
)

reaction(
    label = 'reaction50',
    reactants = ['[CH][O](751)', '[CH]=[C]C=C=[CH](21430)'],
    products = ['[CH]C(C=C=[CH])=C[O](22639)'],
    transitionState = 'TS50',
    kinetics = Arrhenius(A=(2.29708e+06,'m^3/(mol*s)'), n=0.575199, Ea=(34.3157,'kJ/mol'), T0=(1,'K'), comment="""Estimated using an average for rate rule [Y_rad;Birad]
Euclidian distance = 0
Multiplied by reaction path degeneracy 2.0
family: Birad_R_Recombination"""),
)

reaction(
    label = 'reaction51',
    reactants = ['[CH]C(C=C=[CH])=C[O](22639)'],
    products = ['[CH]=[C]C1[CH]C1=C[O](26150)'],
    transitionState = 'TS51',
    kinetics = Arrhenius(A=(2.54e+10,'s^-1'), n=0.69, Ea=(162.626,'kJ/mol'), T0=(1,'K'), comment="""Estimated using an average for rate rule [R4_D_D;doublebond_intra;radadd_intra_cdsingleH]
Euclidian distance = 0
family: Intra_R_Add_Exocyclic
Ea raised from 160.4 to 162.6 kJ/mol to match endothermicity of reaction."""),
)

reaction(
    label = 'reaction52',
    reactants = ['[CH][O](751)', '[CH]=C=CC#C(26151)'],
    products = ['[CH]C(C=C=[CH])=C[O](22639)'],
    transitionState = 'TS52',
    kinetics = Arrhenius(A=(6.03845,'m^3/(mol*s)'), n=1.94267, Ea=(22.8894,'kJ/mol'), T0=(1,'K'), comment="""Estimated using template [Ct-Cd_Ct-H;YJ] for rate rule [Ct-Cd_Ct-H;Y_1centerbirad]
Euclidian distance = 1.0
Multiplied by reaction path degeneracy 2.0
family: R_Addition_MultipleBond"""),
)

reaction(
    label = 'reaction53',
    reactants = ['[CH]=C=C[O](8556)', '[CH]=C=[CH](18734)'],
    products = ['[CH]C(C=C=[CH])=C[O](22639)'],
    transitionState = 'TS53',
    kinetics = Arrhenius(A=(0.246938,'m^3/(mol*s)'), n=2.00579, Ea=(36.0234,'kJ/mol'), T0=(1,'K'), comment="""Estimated using an average for rate rule [Cd_R;CJ]
Euclidian distance = 0
Multiplied by reaction path degeneracy 2.0
family: R_Addition_MultipleBond"""),
)

reaction(
    label = 'reaction54',
    reactants = ['[CH]=C=[C]C([CH2])=C[O](26152)'],
    products = ['[CH]C(C=C=[CH])=C[O](22639)'],
    transitionState = 'TS54',
    kinetics = Arrhenius(A=(4.96e+06,'s^-1'), n=1.85, Ea=(112.257,'kJ/mol'), T0=(1,'K'), comment="""Estimated using template [R3H_SD;Cd_rad_out_double;Cd_H_out_single] for rate rule [R3H_SD;Cd_rad_out_double;Cd_H_out_singleH]
Euclidian distance = 1.0
Multiplied by reaction path degeneracy 2.0
family: intra_H_migration"""),
)

reaction(
    label = 'reaction55',
    reactants = ['[CH]C(C=C=[CH])=C[O](22639)'],
    products = ['[CH]=C1[CH]C([CH]1)=C[O](26092)'],
    transitionState = 'TS55',
    kinetics = Arrhenius(A=(1.953e+11,'s^-1'), n=0.387, Ea=(138.055,'kJ/mol'), T0=(1,'K'), comment="""Estimated using an average for rate rule [R4_D_D;doublebond_intra;radadd_intra_cdsingleH]
Euclidian distance = 0
family: Intra_R_Add_Endocyclic"""),
)

reaction(
    label = 'reaction56',
    reactants = ['[CH][C]=CC([O])C#C(22636)'],
    products = ['[CH]C(C=C=[CH])=C[O](22639)'],
    transitionState = 'TS56',
    kinetics = Arrhenius(A=(3.36138e+09,'s^-1'), n=0.97875, Ea=(163.373,'kJ/mol'), T0=(1,'K'), comment="""Estimated using an average for rate rule [cCsCJ;CJ;C]
Euclidian distance = 0
family: 1,2_shiftC"""),
)

reaction(
    label = 'reaction57',
    reactants = ['[CH]C(C=C=[CH])=C[O](22639)'],
    products = ['C#CC=C1[CH]C1[O](26153)'],
    transitionState = 'TS57',
    kinetics = Arrhenius(A=(2.18e+16,'s^-1'), n=0, Ea=(2.9288,'kJ/mol'), T0=(1,'K'), comment="""Estimated using template [R3;Y_rad_out;Ypri_rad_out] for rate rule [R3_SD;Y_rad_out;CdsinglepriH_rad_out]
Euclidian distance = 2.2360679775
family: Birad_recombination"""),
)

reaction(
    label = 'reaction58',
    reactants = ['H(8)', '[C][C](C=O)C=C=[CH](26154)'],
    products = ['[CH]C(C=C=[CH])=C[O](22639)'],
    transitionState = 'TS58',
    kinetics = Arrhenius(A=(1e+07,'m^3/(mol*s)'), n=0, Ea=(0,'kJ/mol'), T0=(1,'K'), comment="""Estimated using an average for rate rule [H_rad;Birad]
Euclidian distance = 0
family: Birad_R_Recombination"""),
)

reaction(
    label = 'reaction59',
    reactants = ['[CH]C(C=C=[CH])=C[O](22639)'],
    products = ['C#C[CH][C]1[CH]C1[O](26155)'],
    transitionState = 'TS59',
    kinetics = Arrhenius(A=(9.36651e+10,'s^-1'), n=0.5685, Ea=(157.412,'kJ/mol'), T0=(1,'K'), comment="""Estimated using template [R4_D;multiplebond_intra;radadd_intra_cdsingleH] for rate rule [R4_D_CO;carbonylbond_intra_H;radadd_intra_cdsingleH]
Euclidian distance = 2.2360679775
family: Intra_R_Add_Exocyclic
Ea raised from 152.7 to 157.4 kJ/mol to match endothermicity of reaction."""),
)

reaction(
    label = 'reaction60',
    reactants = ['C#CC=O(21959)', '[CH][C]=[CH](21256)'],
    products = ['[CH]C(C=C=[CH])=C[O](22639)'],
    transitionState = 'TS60',
    kinetics = Arrhenius(A=(0.133816,'m^3/(mol*s)'), n=2.39465, Ea=(29.077,'kJ/mol'), T0=(1,'K'), comment="""Estimated using template [Ct-De_Ct-H;CJ] for rate rule [Ct-CO_Ct-H;CJ]
Euclidian distance = 1.0
Multiplied by reaction path degeneracy 2.0
family: R_Addition_MultipleBond"""),
)

reaction(
    label = 'reaction61',
    reactants = ['[CH]=O(373)', '[CH]=[C]C=C=[CH](21430)'],
    products = ['[CH]C(C=C=[CH])=C[O](22639)'],
    transitionState = 'TS61',
    kinetics = Arrhenius(A=(0.133816,'m^3/(mol*s)'), n=2.39465, Ea=(29.077,'kJ/mol'), T0=(1,'K'), comment="""Estimated using template [Ct-De_Ct-H;CJ] for rate rule [Ct-De_Ct-H;CO_pri_rad]
Euclidian distance = 2.0
Multiplied by reaction path degeneracy 2.0
family: R_Addition_MultipleBond"""),
)

reaction(
    label = 'reaction62',
    reactants = ['[CH]C(C=C=[CH])=C[O](22639)'],
    products = ['[CH]=C=C[C]([CH2])[C]=O(26156)'],
    transitionState = 'TS62',
    kinetics = Arrhenius(A=(6718.85,'s^-1'), n=2.58467, Ea=(192.129,'kJ/mol'), T0=(1,'K'), comment="""Estimated using template [R3H_DS;Cd_rad_out_singleH;XH_out] for rate rule [R3H_DS;Cd_rad_out_singleH;CO_H_out]
Euclidian distance = 1.0
family: intra_H_migration"""),
)

reaction(
    label = 'reaction63',
    reactants = ['[CH]=C=C[O](8556)', '[CH][C]=[CH](21256)'],
    products = ['[CH]C(C=C=[CH])=C[O](22639)'],
    transitionState = 'TS63',
    kinetics = Arrhenius(A=(3.9578e+07,'m^3/(mol*s)'), n=-0.126319, Ea=(0,'kJ/mol'), T0=(1,'K'), comment="""Estimated using an average for rate rule [Y_rad;Y_rad]
Euclidian distance = 0
Multiplied by reaction path degeneracy 2.0
family: R_Recombination
Ea raised from -15.6 to -15.6 kJ/mol.
Ea raised from -15.6 to 0 kJ/mol."""),
)

reaction(
    label = 'reaction64',
    reactants = ['[CH]C(C=C=[CH])=C[O](22639)'],
    products = ['[CH]=[C]C=C1[CH]O[CH]1(26157)'],
    transitionState = 'TS64',
    kinetics = Arrhenius(A=(1.953e+11,'s^-1'), n=0.387, Ea=(138.055,'kJ/mol'), T0=(1,'K'), comment="""Estimated using template [R4_D;multiplebond_intra;radadd_intra_cdsingleH] for rate rule [R4_D_CO;carbonyl_intra_H;radadd_intra_cdsingleH]
Euclidian distance = 2.2360679775
family: Intra_R_Add_Endocyclic"""),
)

reaction(
    label = 'reaction65',
    reactants = ['[CH]C(C=C=[CH])=C[O](22639)'],
    products = ['[CH]=C1[CH]C(C=O)=C1(26101)'],
    transitionState = 'TS65',
    kinetics = Arrhenius(A=(1.62e+12,'s^-1'), n=-0.305, Ea=(8.28432,'kJ/mol'), T0=(1,'K'), comment="""Estimated using template [R4;Y_rad_out;Ypri_rad_out] for rate rule [R4_DSD;CdsingleH_rad_out;Ypri_rad_out]
Euclidian distance = 2.2360679775
family: Birad_recombination"""),
)

reaction(
    label = 'reaction21',
    reactants = ['[CH]=[C][CH]OC=C=[CH](22637)'],
    products = ['[CH]C(C=C=[CH])=C[O](22639)'],
    transitionState = 'TS66',
    kinetics = Arrhenius(A=(14080,'s^-1'), n=2.66, Ea=(313.8,'kJ/mol'), T0=(1,'K'), Tmin=(300,'K'), Tmax=(1500,'K'), comment="""Estimated using an average for rate rule [R_ROR;R1_doublebond;R2_doublebond_H;R_O_C]
Euclidian distance = 0
Multiplied by reaction path degeneracy 2.0
family: ketoenol"""),
)

reaction(
    label = 'reaction67',
    reactants = ['[C]#C[CH]C([CH])C=O(26158)'],
    products = ['[CH]C(C=C=[CH])=C[O](22639)'],
    transitionState = 'TS67',
    kinetics = Arrhenius(A=(322699,'s^-1'), n=1.96429, Ea=(152.474,'kJ/mol'), T0=(1,'K'), comment="""Estimated using template [R4Hall;Y_rad_out;XH_out] for rate rule [R4HJ_2;Ct_rad_out;XH_out]
Euclidian distance = 1.41421356237
family: intra_H_migration"""),
)

reaction(
    label = 'reaction68',
    reactants = ['[C]#C[CH]C([CH2])=C[O](26159)'],
    products = ['[CH]C(C=C=[CH])=C[O](22639)'],
    transitionState = 'TS68',
    kinetics = Arrhenius(A=(8820.35,'s^-1'), n=2.27943, Ea=(109.217,'kJ/mol'), T0=(1,'K'), comment="""Estimated using template [R5Hall;Y_rad_out;Cd_H_out_singleH] for rate rule [R5HJ_2;Ct_rad_out;Cd_H_out_singleH]
Euclidian distance = 1.41421356237
Multiplied by reaction path degeneracy 2.0
family: intra_H_migration"""),
)

reaction(
    label = 'reaction69',
    reactants = ['[CH]C(C=C=[CH])=C[O](22639)'],
    products = ['[O]C=C1[CH][C]=C[CH]1(26160)'],
    transitionState = 'TS69',
    kinetics = Arrhenius(A=(2.88e+10,'s^-1'), n=0.31, Ea=(50.6264,'kJ/mol'), T0=(1,'K'), comment="""Estimated using an average for rate rule [R5_linear;triplebond_intra_H;radadd_intra_cdsingleH]
Euclidian distance = 0
family: Intra_R_Add_Endocyclic"""),
)

reaction(
    label = 'reaction70',
    reactants = ['[CH]C(C=C=[CH])=C[O](22639)'],
    products = ['[CH]=C1C([O])C1C#C(26161)'],
    transitionState = 'TS70',
    kinetics = Arrhenius(A=(2e+12,'s^-1'), n=0, Ea=(50.1065,'kJ/mol'), T0=(1,'K'), comment="""Estimated using template [Rn;C_rad_out_H/OneDe;Ypri_rad_out] for rate rule [R3_SS;C_rad_out_H/OneDe;Ypri_rad_out]
Euclidian distance = 2.0
family: Birad_recombination
Ea raised from 48.6 to 50.1 kJ/mol to match endothermicity of reaction."""),
)

reaction(
    label = 'reaction71',
    reactants = ['[CH]C(C=C=[CH])=C[O](22639)'],
    products = ['C#CC1[CH]C1=C[O](26162)'],
    transitionState = 'TS71',
    kinetics = Arrhenius(A=(2e+12,'s^-1'), n=0, Ea=(7.5312,'kJ/mol'), T0=(1,'K'), comment="""Estimated using template [Rn;C_rad_out_H/OneDe;Ypri_rad_out] for rate rule [R3_SD;C_rad_out_H/OneDe;CdsinglepriH_rad_out]
Euclidian distance = 2.82842712475
family: Birad_recombination"""),
)

reaction(
    label = 'reaction72',
    reactants = ['[CH]=[C]C([C]=[CH])C=O(26163)'],
    products = ['[CH]C(C=C=[CH])=C[O](22639)'],
    transitionState = 'TS72',
    kinetics = Arrhenius(A=(1.7778e+12,'s^-1'), n=0.232, Ea=(122.75,'kJ/mol'), T0=(1,'K'), comment="""Estimated using an average for rate rule [cCsCJ;CJ;CO]
Euclidian distance = 0
Multiplied by reaction path degeneracy 2.0
family: 1,2_shiftC"""),
)

network(
    label = '4942',
    isomers = [
        '[CH]C(C=C=[CH])=C[O](22639)',
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
    label = '4942',
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

