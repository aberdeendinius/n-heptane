species(
    label = 'C=[C]C=C[C]=C[O](22614)',
    structure = SMILES('C=[C]C=C[C]=C[O]'),
    E0 = (476.354,'kJ/mol'),
    modes = [
        HarmonicOscillator(frequencies=([2950,3100,1380,975,1025,1650,2995,3010,3025,975,987.5,1000,1300,1337.5,1375,400,450,500,1630,1655,1680,1670,1700,300,440,180,180,180],'cm^-1')),
        HinderedRotor(inertia=(1.63092,'amu*angstrom^2'), symmetry=1, barrier=(37.498,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(1.628,'amu*angstrom^2'), symmetry=1, barrier=(37.4309,'kJ/mol'), semiclassical=False),
    ],
    spinMultiplicity = 4,
    opticalIsomers = 1,
    molecularWeight = (93.1033,'amu'),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[-0.369045,0.0747065,-7.9467e-05,4.04014e-08,-7.54556e-12,57467.9,26.7793], Tmin=(100,'K'), Tmax=(1559.56,'K')), NASAPolynomial(coeffs=[20.0548,0.00664978,1.06533e-06,-4.6779e-10,3.88295e-14,53003.4,-74.7195], Tmin=(1559.56,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(476.354,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(274.378,'J/(mol*K)'), comment="""Thermo group additivity estimation: group(O2s-(Cds-Cd)H) + group(Cds-Cds(Cds-Cds)H) + group(Cds-Cds(Cds-Cds)H) + group(Cds-Cds(Cds-Cds)H) + group(Cds-Cds(Cds-Cds)H) + group(Cds-CdsOsH) + group(Cds-CdsHH) + radical(C=COJ) + radical(C=CJC=C) + radical(C=CJC=C)"""),
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
    label = 'C=[C][CH]C1[C]=CO1(24975)',
    structure = SMILES('C=[C][CH]C1[C]=CO1'),
    E0 = (588.552,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
    molecularWeight = (93.1033,'amu'),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[1.39124,0.0271944,8.17993e-05,-1.45444e-07,6.28578e-11,70908.3,23.8494], Tmin=(100,'K'), Tmax=(930.985,'K')), NASAPolynomial(coeffs=[26.279,-0.00231864,4.61514e-06,-8.51852e-10,4.43407e-14,62919.2,-112.459], Tmin=(930.985,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(588.552,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(274.378,'J/(mol*K)'), comment="""Thermo group additivity estimation: group(O2s-Cs(Cds-Cd)) + group(Cs-(Cds-Cds)CsOsH) + group(Cs-(Cds-Cds)CsHH) + group(Cds-CdsCsH) + group(Cds-CdsCsH) + group(Cds-CdsOsH) + group(Cds-CdsHH) + ring(Cyclobutene) + radical(Cds_S) + radical(Cds_S) + radical(C=CCJC(O)C=C)"""),
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
    label = '[CH2][C]=CC#CC=O(24977)',
    structure = SMILES('[CH2][C]=CC#CC=O'),
    E0 = (463.979,'kJ/mol'),
    modes = [
        HarmonicOscillator(frequencies=([1685,370,2782.5,750,1395,475,1775,1000,3010,987.5,1337.5,450,1655,2100,2250,500,550,3000,3100,440,815,1455,1000,236.731],'cm^-1')),
        HinderedRotor(inertia=(0.515097,'amu*angstrom^2'), symmetry=1, barrier=(20.5018,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.979099,'amu*angstrom^2'), symmetry=1, barrier=(38.933,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(1.80648,'amu*angstrom^2'), symmetry=1, barrier=(71.8069,'kJ/mol'), semiclassical=False),
    ],
    spinMultiplicity = 3,
    opticalIsomers = 1,
    molecularWeight = (92.0954,'amu'),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[1.82884,0.0498187,-4.53977e-05,2.29638e-08,-4.92108e-12,55880.3,21.445], Tmin=(100,'K'), Tmax=(1089.09,'K')), NASAPolynomial(coeffs=[8.53476,0.0251894,-1.14761e-05,2.19956e-09,-1.54694e-13,54419.6,-11.4796], Tmin=(1089.09,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(463.979,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(245.277,'J/(mol*K)'), comment="""Thermo group additivity estimation: group(Cs-(Cds-Cds)HHH) + group(Cds-CdsCsH) + group(Cds-CdsCtH) + group(Cds-O2d(Cds-Cds)H) + group(Ct-Ct(Cds-Cds)) + group(Ct-CtCs) + radical(CTCC=CCJ) + radical(Cds_S)"""),
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
    label = 'C=C=C[CH][C]=C=O(24978)',
    structure = SMILES('C=C=C[CH][C]=C=O'),
    E0 = (446.256,'kJ/mol'),
    modes = [
        HarmonicOscillator(frequencies=([3025,407.5,1350,352.5,540,610,2055,2950,3100,1380,975,1025,1650,2120,512.5,787.5,1685,370,3010,987.5,1337.5,450,1655,180,180],'cm^-1')),
        HinderedRotor(inertia=(1.7166,'amu*angstrom^2'), symmetry=1, barrier=(39.4679,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(1.714,'amu*angstrom^2'), symmetry=1, barrier=(39.4081,'kJ/mol'), semiclassical=False),
    ],
    spinMultiplicity = 3,
    opticalIsomers = 1,
    molecularWeight = (92.0954,'amu'),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[1.19851,0.0571887,-5.44572e-05,2.64267e-08,-5.09381e-12,53776.8,21.2307], Tmin=(100,'K'), Tmax=(1254.7,'K')), NASAPolynomial(coeffs=[13.5989,0.0176565,-7.19652e-06,1.31559e-09,-9.04376e-14,50665,-41.4076], Tmin=(1254.7,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(446.256,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(249.434,'J/(mol*K)'), comment="""Thermo group additivity estimation: group(Cs-Cd(CCO)HH) + group(Cds-CdsCsH) + group(Cds-(Cdd-O2d)CsH) + group(Cds-CdsHH) + group(Cdd-CdsCds) + radical(CCCJ=C=O) + radical(C=CCJC=C=O)"""),
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
    label = 'C#CC=C[C]=C(21245)',
    structure = SMILES('C#CC=C[C]=C'),
    E0 = (521.475,'kJ/mol'),
    modes = [
        HarmonicOscillator(frequencies=([2950,3100,1380,975,1025,1650,750,770,3400,2100,2995,3025,975,1000,1300,1375,400,500,1630,1680,2175,525,1685,370,180],'cm^-1')),
        HinderedRotor(inertia=(1.62249,'amu*angstrom^2'), symmetry=1, barrier=(37.3042,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(1.61979,'amu*angstrom^2'), symmetry=1, barrier=(37.2421,'kJ/mol'), semiclassical=False),
    ],
    spinMultiplicity = 2,
    opticalIsomers = 1,
    molecularWeight = (77.1039,'amu'),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[1.45098,0.0475986,-2.51853e-05,-8.16809e-09,8.69319e-12,62818.4,18.2704], Tmin=(100,'K'), Tmax=(937.719,'K')), NASAPolynomial(coeffs=[14.4661,0.0131635,-3.8272e-06,6.24214e-10,-4.31862e-14,59450.6,-48.6255], Tmin=(937.719,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(521.475,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(249.434,'J/(mol*K)'), comment="""Thermo group additivity estimation: group(Cds-Cds(Cds-Cds)H) + group(Cds-Cds(Cds-Cds)H) + group(Cds-CdsCtH) + group(Cds-CdsHH) + group(Ct-Ct(Cds-Cds)) + group(Ct-CtH) + radical(C=CJC=C)"""),
)

species(
    label = '[CH2]C=[C]C=C=C[O](24979)',
    structure = SMILES('[CH2]C=[C]C=C=C[O]'),
    E0 = (448.196,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
    molecularWeight = (93.1033,'amu'),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[0.987826,0.0543355,-2.02451e-05,-2.49503e-08,1.75151e-11,54025.1,23.8978], Tmin=(100,'K'), Tmax=(909.932,'K')), NASAPolynomial(coeffs=[18.1116,0.011083,-1.73239e-06,1.61599e-10,-1.02127e-14,49583.1,-64.3826], Tmin=(909.932,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(448.196,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(274.378,'J/(mol*K)'), comment="""Thermo group additivity estimation: group(O2s-(Cds-Cd)H) + group(Cs-(Cds-Cds)HHH) + group(Cds-CdsCsH) + group(Cds-Cds(Cds-Cds)H) + group(Cds-Cds(Cds-Cds)H) + group(Cds-CdsOsH) + group(Cdd-CdsCds) + radical(C=CJC=C) + radical(C=CC=CCJ) + radical(C=COJ)"""),
)

species(
    label = 'C=[C]C=[C]C=C[O](24980)',
    structure = SMILES('C=[C]C=[C]C=C[O]'),
    E0 = (476.354,'kJ/mol'),
    modes = [
        HarmonicOscillator(frequencies=([2950,3100,1380,975,1025,1650,2995,3010,3025,975,987.5,1000,1300,1337.5,1375,400,450,500,1630,1655,1680,1670,1700,300,440,180,180,180],'cm^-1')),
        HinderedRotor(inertia=(1.63092,'amu*angstrom^2'), symmetry=1, barrier=(37.498,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(1.628,'amu*angstrom^2'), symmetry=1, barrier=(37.4309,'kJ/mol'), semiclassical=False),
    ],
    spinMultiplicity = 4,
    opticalIsomers = 1,
    molecularWeight = (93.1033,'amu'),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[-0.369045,0.0747065,-7.9467e-05,4.04014e-08,-7.54556e-12,57467.9,26.7793], Tmin=(100,'K'), Tmax=(1559.56,'K')), NASAPolynomial(coeffs=[20.0548,0.00664978,1.06533e-06,-4.6779e-10,3.88295e-14,53003.4,-74.7195], Tmin=(1559.56,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(476.354,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(274.378,'J/(mol*K)'), comment="""Thermo group additivity estimation: group(O2s-(Cds-Cd)H) + group(Cds-Cds(Cds-Cds)H) + group(Cds-Cds(Cds-Cds)H) + group(Cds-Cds(Cds-Cds)H) + group(Cds-Cds(Cds-Cds)H) + group(Cds-CdsOsH) + group(Cds-CdsHH) + radical(C=CJC=C) + radical(C=COJ) + radical(C=CJC=C)"""),
)

species(
    label = '[CH]C=CC=C=C[O](24787)',
    structure = SMILES('[CH]C=CC=C=C[O]'),
    E0 = (501.829,'kJ/mol'),
    modes = [
        HarmonicOscillator(frequencies=([2995,3005,3015,3025,975,983.333,991.667,1000,1300,1325,1350,1375,400,433.333,466.667,500,1630,1646.67,1663.33,1680,540,610,2055,180,180,180,180,180],'cm^-1')),
        HinderedRotor(inertia=(1.99804,'amu*angstrom^2'), symmetry=1, barrier=(45.9389,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(1.99879,'amu*angstrom^2'), symmetry=1, barrier=(45.956,'kJ/mol'), semiclassical=False),
    ],
    spinMultiplicity = 4,
    opticalIsomers = 1,
    molecularWeight = (93.1033,'amu'),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[0.782723,0.0566767,-1.4346e-05,-3.18047e-08,1.91222e-11,60484.8,24.5867], Tmin=(100,'K'), Tmax=(940.73,'K')), NASAPolynomial(coeffs=[18.7528,0.0153334,-4.3366e-06,7.25523e-10,-5.27467e-14,55552.2,-69.2575], Tmin=(940.73,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(501.829,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(274.378,'J/(mol*K)'), comment="""Thermo group additivity estimation: group(O2s-(Cds-Cd)H) + group(Cs-(Cds-Cds)HHH) + group(Cds-CdsCsH) + group(Cds-Cds(Cds-Cds)H) + group(Cds-Cds(Cds-Cds)H) + group(Cds-CdsOsH) + group(Cdd-CdsCds) + radical(AllylJ2_triplet) + radical(C=COJ)"""),
)

species(
    label = 'C=[C]C=C[C]=[C]O(24981)',
    structure = SMILES('C=[C]C=C[C]=[C]O'),
    E0 = (574.636,'kJ/mol'),
    modes = [
        HarmonicOscillator(frequencies=([1670,1685,1700,300,370,440,2950,3100,1380,975,1025,1650,3615,1277.5,1000,2995,3025,975,1000,1300,1375,400,500,1630,1680,180,180],'cm^-1')),
        HinderedRotor(inertia=(1.39948,'amu*angstrom^2'), symmetry=1, barrier=(32.1769,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(1.39571,'amu*angstrom^2'), symmetry=1, barrier=(32.0902,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(1.39965,'amu*angstrom^2'), symmetry=1, barrier=(32.1806,'kJ/mol'), semiclassical=False),
    ],
    spinMultiplicity = 4,
    opticalIsomers = 1,
    molecularWeight = (93.1033,'amu'),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[-0.144368,0.0786389,-9.39946e-05,5.38637e-08,-1.14118e-11,69272.8,27.6757], Tmin=(100,'K'), Tmax=(1344.25,'K')), NASAPolynomial(coeffs=[19.672,0.00658976,1.00128e-06,-4.88604e-10,4.30368e-14,65127.2,-69.3924], Tmin=(1344.25,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(574.636,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(270.22,'J/(mol*K)'), comment="""Thermo group additivity estimation: group(O2s-(Cds-Cd)H) + group(Cds-Cds(Cds-Cds)H) + group(Cds-Cds(Cds-Cds)H) + group(Cds-Cds(Cds-Cds)H) + group(Cds-Cds(Cds-Cds)H) + group(Cds-CdsOsH) + group(Cds-CdsHH) + radical(C=CJC=C) + radical(C=CJO) + radical(C=CJC=C)"""),
)

species(
    label = '[CH2][C]=CC=C[C]=O(24982)',
    structure = SMILES('[CH2][C]=CC=C[C]=O'),
    E0 = (450.927,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
    molecularWeight = (93.1033,'amu'),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[1.25041,0.0629723,-6.63244e-05,3.72924e-08,-8.59495e-12,54330.9,24.3127], Tmin=(100,'K'), Tmax=(1037,'K')), NASAPolynomial(coeffs=[10.9988,0.0253702,-1.1934e-05,2.32611e-09,-1.65302e-13,52309.1,-23.0719], Tmin=(1037,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(450.927,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(270.22,'J/(mol*K)'), comment="""Thermo group additivity estimation: group(Cs-(Cds-Cds)HHH) + group(Cds-CdsCsH) + group(Cds-Cds(Cds-Cds)H) + group(Cds-Cds(Cds-Cds)H) + group(Cd-Cd(CO)H) + group(Cds-O2d(Cds-Cds)H) + radical(C=CC=CCJ) + radical(Cds_S) + radical(C=CCJ=O)"""),
)

species(
    label = 'C=CC=[C][C]=C[O](24983)',
    structure = SMILES('C=CC=[C][C]=C[O]'),
    E0 = (476.354,'kJ/mol'),
    modes = [
        HarmonicOscillator(frequencies=([2950,3100,1380,975,1025,1650,2995,3010,3025,975,987.5,1000,1300,1337.5,1375,400,450,500,1630,1655,1680,1670,1700,300,440,180,180,180],'cm^-1')),
        HinderedRotor(inertia=(1.63092,'amu*angstrom^2'), symmetry=1, barrier=(37.498,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(1.628,'amu*angstrom^2'), symmetry=1, barrier=(37.4309,'kJ/mol'), semiclassical=False),
    ],
    spinMultiplicity = 4,
    opticalIsomers = 1,
    molecularWeight = (93.1033,'amu'),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[-0.369045,0.0747065,-7.9467e-05,4.04014e-08,-7.54556e-12,57467.9,26.7793], Tmin=(100,'K'), Tmax=(1559.56,'K')), NASAPolynomial(coeffs=[20.0548,0.00664978,1.06533e-06,-4.6779e-10,3.88295e-14,53003.4,-74.7195], Tmin=(1559.56,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(476.354,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(274.378,'J/(mol*K)'), comment="""Thermo group additivity estimation: group(O2s-(Cds-Cd)H) + group(Cds-Cds(Cds-Cds)H) + group(Cds-Cds(Cds-Cds)H) + group(Cds-Cds(Cds-Cds)H) + group(Cds-Cds(Cds-Cds)H) + group(Cds-CdsOsH) + group(Cds-CdsHH) + radical(C=COJ) + radical(C=CJC=C) + radical(C=CJC=C)"""),
)

species(
    label = 'C=[C][C]=CC=C[O](24984)',
    structure = SMILES('C=[C][C]=CC=C[O]'),
    E0 = (476.354,'kJ/mol'),
    modes = [
        HarmonicOscillator(frequencies=([2950,3100,1380,975,1025,1650,2995,3010,3025,975,987.5,1000,1300,1337.5,1375,400,450,500,1630,1655,1680,1670,1700,300,440,180,180,180],'cm^-1')),
        HinderedRotor(inertia=(1.63092,'amu*angstrom^2'), symmetry=1, barrier=(37.498,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(1.628,'amu*angstrom^2'), symmetry=1, barrier=(37.4309,'kJ/mol'), semiclassical=False),
    ],
    spinMultiplicity = 4,
    opticalIsomers = 1,
    molecularWeight = (93.1033,'amu'),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[-0.369045,0.0747065,-7.9467e-05,4.04014e-08,-7.54556e-12,57467.9,26.7793], Tmin=(100,'K'), Tmax=(1559.56,'K')), NASAPolynomial(coeffs=[20.0548,0.00664978,1.06533e-06,-4.6779e-10,3.88295e-14,53003.4,-74.7195], Tmin=(1559.56,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(476.354,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(274.378,'J/(mol*K)'), comment="""Thermo group additivity estimation: group(O2s-(Cds-Cd)H) + group(Cds-Cds(Cds-Cds)H) + group(Cds-Cds(Cds-Cds)H) + group(Cds-Cds(Cds-Cds)H) + group(Cds-Cds(Cds-Cds)H) + group(Cds-CdsOsH) + group(Cds-CdsHH) + radical(C=COJ) + radical(C=CJC=C) + radical(C=CJC=C)"""),
)

species(
    label = 'C=[C]C=[C][C]=CO(24985)',
    structure = SMILES('C=[C]C=[C][C]=CO'),
    E0 = (533.887,'kJ/mol'),
    modes = [
        HarmonicOscillator(frequencies=([1670,1685,1700,300,370,440,2950,3100,1380,975,1025,1650,3615,1277.5,1000,2995,3025,975,1000,1300,1375,400,500,1630,1680,180,180],'cm^-1')),
        HinderedRotor(inertia=(1.61181,'amu*angstrom^2'), symmetry=1, barrier=(37.0588,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(1.60975,'amu*angstrom^2'), symmetry=1, barrier=(37.0114,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(1.61996,'amu*angstrom^2'), symmetry=1, barrier=(37.2461,'kJ/mol'), semiclassical=False),
    ],
    spinMultiplicity = 4,
    opticalIsomers = 1,
    molecularWeight = (93.1033,'amu'),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[-0.6451,0.086018,-0.000105693,6.09116e-08,-1.27855e-11,64393,25.9355], Tmin=(100,'K'), Tmax=(1389.15,'K')), NASAPolynomial(coeffs=[21.5834,0.0035262,3.34205e-06,-9.95048e-10,7.95649e-14,60000.9,-82.1904], Tmin=(1389.15,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(533.887,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(270.22,'J/(mol*K)'), comment="""Thermo group additivity estimation: group(O2s-(Cds-Cd)H) + group(Cds-Cds(Cds-Cds)H) + group(Cds-Cds(Cds-Cds)H) + group(Cds-Cds(Cds-Cds)H) + group(Cds-Cds(Cds-Cds)H) + group(Cds-CdsOsH) + group(Cds-CdsHH) + radical(C=CJC=C) + radical(C=CJC=C) + radical(C=CJC=C)"""),
)

species(
    label = '[CH]=C=C[CH]C=C[O](21998)',
    structure = SMILES('[CH]=C=C[CH]C=C[O]'),
    E0 = (418.2,'kJ/mol'),
    modes = [
        HarmonicOscillator(frequencies=([3120,650,792.5,1650,3025,407.5,1350,352.5,2995,3010,3025,975,987.5,1000,1300,1337.5,1375,400,450,500,1630,1655,1680,540,610,2055,180,180],'cm^-1')),
        HinderedRotor(inertia=(1.81308,'amu*angstrom^2'), symmetry=1, barrier=(41.6863,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(1.81229,'amu*angstrom^2'), symmetry=1, barrier=(41.668,'kJ/mol'), semiclassical=False),
    ],
    spinMultiplicity = 4,
    opticalIsomers = 1,
    molecularWeight = (93.1033,'amu'),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[1.33213,0.0407037,2.28675e-05,-7.3903e-08,3.61183e-11,50410.7,24.8266], Tmin=(100,'K'), Tmax=(908.639,'K')), NASAPolynomial(coeffs=[20.0131,0.00610111,1.35435e-06,-4.23921e-10,2.75502e-14,45049.5,-74.3289], Tmin=(908.639,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(418.2,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(274.378,'J/(mol*K)'), comment="""Thermo group additivity estimation: group(O2s-(Cds-Cd)H) + group(Cs-(Cds-Cds)(Cds-Cds)HH) + group(Cds-CdsCsH) + group(Cds-CdsCsH) + group(Cds-CdsOsH) + group(Cds-CdsHH) + group(Cdd-CdsCds) + radical(C=CCJC=C) + radical(C=COJ) + radical(C=C=CJ)"""),
)

species(
    label = '[CH2]C=C[CH][C]=C=O(24986)',
    structure = SMILES('[CH2]C=C[CH][C]=C=O'),
    E0 = (421.151,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
    molecularWeight = (93.1033,'amu'),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[1.0914,0.0558621,-3.97911e-05,8.62889e-09,1.43592e-12,50764.5,22.8412], Tmin=(100,'K'), Tmax=(1035.59,'K')), NASAPolynomial(coeffs=[14.4684,0.018809,-7.292e-06,1.33617e-09,-9.36624e-14,47210.2,-45.9469], Tmin=(1035.59,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(421.151,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(270.22,'J/(mol*K)'), comment="""Thermo group additivity estimation: group(Cs-Cd(CCO)HH) + group(Cs-(Cds-Cds)HHH) + group(Cds-CdsCsH) + group(Cds-CdsCsH) + group(Cds-(Cdd-O2d)CsH) + radical(CCCJ=C=O) + radical(C=CCJC=C=O) + radical(Allyl_P)"""),
)

species(
    label = 'C=[C][C]=C[C]=CO(24987)',
    structure = SMILES('C=[C][C]=C[C]=CO'),
    E0 = (533.887,'kJ/mol'),
    modes = [
        HarmonicOscillator(frequencies=([1670,1685,1700,300,370,440,2950,3100,1380,975,1025,1650,3615,1277.5,1000,2995,3025,975,1000,1300,1375,400,500,1630,1680,180,180],'cm^-1')),
        HinderedRotor(inertia=(1.61181,'amu*angstrom^2'), symmetry=1, barrier=(37.0588,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(1.60975,'amu*angstrom^2'), symmetry=1, barrier=(37.0114,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(1.61996,'amu*angstrom^2'), symmetry=1, barrier=(37.2461,'kJ/mol'), semiclassical=False),
    ],
    spinMultiplicity = 4,
    opticalIsomers = 1,
    molecularWeight = (93.1033,'amu'),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[-0.6451,0.086018,-0.000105693,6.09116e-08,-1.27855e-11,64393,25.9355], Tmin=(100,'K'), Tmax=(1389.15,'K')), NASAPolynomial(coeffs=[21.5834,0.0035262,3.34205e-06,-9.95048e-10,7.95649e-14,60000.9,-82.1904], Tmin=(1389.15,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(533.887,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(270.22,'J/(mol*K)'), comment="""Thermo group additivity estimation: group(O2s-(Cds-Cd)H) + group(Cds-Cds(Cds-Cds)H) + group(Cds-Cds(Cds-Cds)H) + group(Cds-Cds(Cds-Cds)H) + group(Cds-Cds(Cds-Cds)H) + group(Cds-CdsOsH) + group(Cds-CdsHH) + radical(C=CJC=C) + radical(C=CJC=C) + radical(C=CJC=C)"""),
)

species(
    label = '[CH]=C=C[CH][C]=CO(24792)',
    structure = SMILES('[CH]=C=C[CH][C]=CO'),
    E0 = (514.579,'kJ/mol'),
    modes = [
        HarmonicOscillator(frequencies=([3120,650,792.5,1650,3615,1277.5,1000,3025,407.5,1350,352.5,540,610,2055,1685,370,2995,3025,975,1000,1300,1375,400,500,1630,1680,180],'cm^-1')),
        HinderedRotor(inertia=(1.45877,'amu*angstrom^2'), symmetry=1, barrier=(33.5399,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(1.45527,'amu*angstrom^2'), symmetry=1, barrier=(33.4595,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(1.45433,'amu*angstrom^2'), symmetry=1, barrier=(33.4379,'kJ/mol'), semiclassical=False),
    ],
    spinMultiplicity = 4,
    opticalIsomers = 1,
    molecularWeight = (93.1033,'amu'),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[0.892675,0.0521838,-4.83939e-06,-5.17907e-08,3.04716e-11,62016.8,25.4694], Tmin=(100,'K'), Tmax=(889.314,'K')), NASAPolynomial(coeffs=[22.1894,0.00172589,3.80739e-06,-9.54788e-10,6.78885e-14,56436.3,-84.8552], Tmin=(889.314,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(514.579,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(270.22,'J/(mol*K)'), comment="""Thermo group additivity estimation: group(O2s-(Cds-Cd)H) + group(Cs-(Cds-Cds)(Cds-Cds)HH) + group(Cds-CdsCsH) + group(Cds-CdsCsH) + group(Cds-CdsOsH) + group(Cds-CdsHH) + group(Cdd-CdsCds) + radical(Cds_S) + radical(C=CCJC=C) + radical(C=C=CJ)"""),
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
    label = 'C=C1[CH]C1[C]=C[O](24883)',
    structure = SMILES('C=C1[CH]C1[C]=C[O]'),
    E0 = (570.271,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
    molecularWeight = (93.1033,'amu'),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[1.25758,0.047576,-8.87282e-06,-3.01777e-08,1.70635e-11,68698,20.2945], Tmin=(100,'K'), Tmax=(952.453,'K')), NASAPolynomial(coeffs=[16.8363,0.0132953,-3.93459e-06,6.98075e-10,-5.22834e-14,64317.7,-61.5211], Tmin=(952.453,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(570.271,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(278.535,'J/(mol*K)'), comment="""Thermo group additivity estimation: group(O2s-(Cds-Cd)H) + group(Cs-(Cds-Cds)(Cds-Cds)CsH) + group(Cs-(Cds-Cds)CsHH) + group(Cds-CdsCsCs) + group(Cds-CdsCsH) + group(Cds-CdsOsH) + group(Cds-CdsHH) + ring(Methylene_cyclopropane) + radical(Cds_S) + radical(C=COJ) + radical(Allyl_S)"""),
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
    label = 'C=C=C[CH]C=C=O(24992)',
    structure = SMILES('C=C=C[CH]C=C=O'),
    E0 = (243.975,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
    molecularWeight = (93.1033,'amu'),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[1.15368,0.0568041,-4.92735e-05,2.18828e-08,-3.88598e-12,29450.7,21.6849], Tmin=(100,'K'), Tmax=(1350.35,'K')), NASAPolynomial(coeffs=[13.7544,0.0194783,-7.81116e-06,1.41285e-09,-9.62338e-14,26047.7,-42.8911], Tmin=(1350.35,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(243.975,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(274.378,'J/(mol*K)'), comment="""Thermo group additivity estimation: group(Cs-Cd(CCO)HH) + group(Cds-CdsCsH) + group(Cds-(Cdd-O2d)CsH) + group(Cds-CdsHH) + group(Cdd-CdsCds) + radical(C=CCJC=C=O)"""),
)

species(
    label = '[CH2]C#CC=CC=O(24993)',
    structure = SMILES('[CH2]C#CC=CC=O'),
    E0 = (261.445,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
    molecularWeight = (93.1033,'amu'),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[1.55419,0.0499545,-3.48774e-05,1.12288e-08,-1.4285e-12,31535.4,23.2374], Tmin=(100,'K'), Tmax=(1812.95,'K')), NASAPolynomial(coeffs=[15.6984,0.0187474,-9.05718e-06,1.73399e-09,-1.19198e-13,26406.9,-53.4153], Tmin=(1812.95,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(261.445,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(270.22,'J/(mol*K)'), comment="""Thermo group additivity estimation: group(Cs-CtHHH) + group(Cd-Cd(CO)H) + group(Cds-CdsCtH) + group(Cds-O2d(Cds-Cds)H) + group(Ct-CtCs) + group(Ct-Ct(Cds-Cds)) + radical(Propargyl)"""),
)

species(
    label = '[CH2]C=CC#CC=O(24994)',
    structure = SMILES('[CH2]C=CC#CC=O'),
    E0 = (226.137,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
    molecularWeight = (93.1033,'amu'),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[1.45307,0.049942,-3.65993e-05,1.3567e-08,-2.03487e-12,27294.7,22.3247], Tmin=(100,'K'), Tmax=(1558.94,'K')), NASAPolynomial(coeffs=[12.8707,0.0206456,-8.40991e-06,1.51183e-09,-1.01602e-13,23734.9,-37.828], Tmin=(1558.94,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(226.137,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(270.22,'J/(mol*K)'), comment="""Thermo group additivity estimation: group(Cs-(Cds-Cds)HHH) + group(Cds-CdsCsH) + group(Cds-CdsCtH) + group(Cds-O2d(Cds-Cds)H) + group(Ct-Ct(Cds-Cds)) + group(Ct-CtCs) + radical(CTCC=CCJ)"""),
)

species(
    label = 'C=[C]C=C=C=CO(24995)',
    structure = SMILES('C=[C]C=C=C=CO'),
    E0 = (365.281,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
    molecularWeight = (93.1033,'amu'),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[-0.470013,0.0793851,-8.95316e-05,4.72805e-08,-9.17474e-12,44110.4,25.2131], Tmin=(100,'K'), Tmax=(1478.64,'K')), NASAPolynomial(coeffs=[21.9953,0.00400058,1.76463e-06,-5.64834e-10,4.46129e-14,39064.1,-86.5546], Tmin=(1478.64,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(365.281,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(274.378,'J/(mol*K)'), comment="""Thermo group additivity estimation: group(O2s-(Cds-Cd)H) + group(Cds-Cds(Cds-Cds)H) + group(Cds-Cds(Cds-Cds)H) + group(Cds-CdsHH) + group(Cds-CdsOsH) + group(Cdd-CdsCds) + group(Cdd-CdsCds) + radical(C=CJC=C)"""),
)

species(
    label = 'C=C=CC=C1[CH]O1(24996)',
    structure = SMILES('C=C=CC=C1[CH]O1'),
    E0 = (333.931,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
    molecularWeight = (93.1033,'amu'),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[1.20971,0.0355861,5.21806e-05,-1.1374e-07,5.19108e-11,40287,19.5808], Tmin=(100,'K'), Tmax=(922.364,'K')), NASAPolynomial(coeffs=[25.3042,-0.00209112,4.79883e-06,-9.60115e-10,5.66978e-14,33000.2,-110.12], Tmin=(922.364,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(333.931,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(278.535,'J/(mol*K)'), comment="""Thermo group additivity estimation: group(O2s-Cs(Cds-Cd)) + group(Cs-(Cds-Cds)OsHH) + group(Cds-CdsCsOs) + group(Cds-Cds(Cds-Cds)H) + group(Cds-Cds(Cds-Cds)H) + group(Cds-CdsHH) + group(Cdd-CdsCds) + ring(methyleneoxirane) + radical(C=CCJO)"""),
)

species(
    label = 'C=C1C=C[C]1C=O(24864)',
    structure = SMILES('C=C1C=C[C]1C=O'),
    E0 = (195.27,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
    molecularWeight = (93.1033,'amu'),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[1.92718,0.0269745,4.96814e-05,-8.82455e-08,3.67333e-11,23577.1,16.7175], Tmin=(100,'K'), Tmax=(953.721,'K')), NASAPolynomial(coeffs=[15.9961,0.0146426,-4.33234e-06,8.25221e-10,-6.59894e-14,18770.8,-61.6191], Tmin=(953.721,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(195.27,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(278.535,'J/(mol*K)'), comment="""Thermo group additivity estimation: group(Cs-(Cds-O2d)CsCsH) + group(Cds-Cds(Cds-Cds)Cs) + group(Cds-CdsCsH) + group(Cds-Cds(Cds-Cds)H) + group(Cds-OdCsH) + group(Cds-CdsHH) + ring(Cyclobutene) + radical(C=CCJ(C=O)C=C)"""),
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
    label = '[CH]C=C=C(5316)',
    structure = SMILES('[CH]C=C=C'),
    E0 = (517.387,'kJ/mol'),
    modes = [
        HarmonicOscillator(frequencies=([2950,3100,1380,975,1025,1650,3010,987.5,1337.5,450,1655,540,610,2055,180,180,180],'cm^-1')),
        HinderedRotor(inertia=(2.12858,'amu*angstrom^2'), symmetry=1, barrier=(48.9403,'kJ/mol'), semiclassical=False),
    ],
    spinMultiplicity = 3,
    opticalIsomers = 1,
    molecularWeight = (52.0746,'amu'),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[2.68278,0.0236277,3.95557e-06,-1.71396e-08,6.98364e-12,62279.3,13.8634], Tmin=(100,'K'), Tmax=(1041.69,'K')), NASAPolynomial(coeffs=[6.67252,0.0201465,-8.07948e-06,1.47317e-09,-1.01811e-13,60805.7,-8.63085], Tmin=(1041.69,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(517.387,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(178.761,'J/(mol*K)'), comment="""Thermo group additivity estimation: group(Cs-(Cds-Cds)HHH) + group(Cds-CdsCsH) + group(Cds-CdsHH) + group(Cdd-CdsCds) + radical(AllylJ2_triplet)"""),
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
    label = '[CH]=C=CC=[C][CH2](19899)',
    structure = SMILES('[CH]=C=CC=[C][CH2]'),
    E0 = (708.849,'kJ/mol'),
    modes = [
        HarmonicOscillator(frequencies=([3120,650,792.5,1650,1685,370,540,610,2055,2995,3025,975,1000,1300,1375,400,500,1630,1680,3000,3100,440,815,1455,1000],'cm^-1')),
        HinderedRotor(inertia=(1.83532,'amu*angstrom^2'), symmetry=1, barrier=(42.1976,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(1.83581,'amu*angstrom^2'), symmetry=1, barrier=(42.2088,'kJ/mol'), semiclassical=False),
    ],
    spinMultiplicity = 4,
    opticalIsomers = 1,
    molecularWeight = (77.1039,'amu'),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[1.11113,0.051108,-4.51432e-05,2.08072e-08,-3.69078e-12,85369.4,23.6956], Tmin=(100,'K'), Tmax=(1566.08,'K')), NASAPolynomial(coeffs=[13.2742,0.0137115,-3.26158e-06,3.97592e-10,-2.06981e-14,82336,-37.962], Tmin=(1566.08,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(708.849,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(249.434,'J/(mol*K)'), comment="""Thermo group additivity estimation: group(Cs-(Cds-Cds)HHH) + group(Cds-CdsCsH) + group(Cds-Cds(Cds-Cds)H) + group(Cds-Cds(Cds-Cds)H) + group(Cds-CdsHH) + group(Cdd-CdsCds) + radical(C=C=CJ) + radical(C=CC=CCJ) + radical(Cds_S)"""),
)

species(
    label = '[O]C=[C]C1C=[C]C1(24785)',
    structure = SMILES('[O]C=[C]C1C=[C]C1'),
    E0 = (639.167,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
    molecularWeight = (93.1033,'amu'),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[1.36876,0.0465842,-1.22713e-05,-2.28587e-08,1.3543e-11,76978.9,22.5998], Tmin=(100,'K'), Tmax=(965.255,'K')), NASAPolynomial(coeffs=[15.4347,0.0151074,-5.02228e-06,9.11612e-10,-6.66601e-14,73014.4,-51.2328], Tmin=(965.255,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(639.167,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(278.535,'J/(mol*K)'), comment="""Thermo group additivity estimation: group(O2s-(Cds-Cd)H) + group(Cs-(Cds-Cds)(Cds-Cds)CsH) + group(Cs-(Cds-Cds)CsHH) + group(Cds-CdsCsH) + group(Cds-CdsCsH) + group(Cds-CdsCsH) + group(Cds-CdsOsH) + ring(Cyclobutene) + radical(Cds_S) + radical(C=COJ) + radical(cyclobutene-vinyl)"""),
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
    label = 'C[C]=[C]C=C=C[O](24999)',
    structure = SMILES('C[C]=[C]C=C=C[O]'),
    E0 = (567.982,'kJ/mol'),
    modes = [
        HarmonicOscillator(frequencies=([2750,2800,2850,1350,1500,750,1050,1375,1000,540,610,2055,2995,3025,975,1000,1300,1375,400,500,1630,1680,1670,1700,300,440,180,180],'cm^-1')),
        HinderedRotor(inertia=(1.36921,'amu*angstrom^2'), symmetry=1, barrier=(31.4809,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(1.37912,'amu*angstrom^2'), symmetry=1, barrier=(31.7087,'kJ/mol'), semiclassical=False),
    ],
    spinMultiplicity = 4,
    opticalIsomers = 1,
    molecularWeight = (93.1033,'amu'),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[0.451466,0.0705229,-7.79491e-05,4.32567e-08,-9.24688e-12,68446.6,24.7351], Tmin=(100,'K'), Tmax=(1195.86,'K')), NASAPolynomial(coeffs=[16.9644,0.0133388,-3.7752e-06,5.42506e-10,-3.21616e-14,64636.6,-57.3008], Tmin=(1195.86,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(567.982,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(274.378,'J/(mol*K)'), comment="""Thermo group additivity estimation: group(O2s-(Cds-Cd)H) + group(Cs-(Cds-Cds)HHH) + group(Cds-CdsCsH) + group(Cds-Cds(Cds-Cds)H) + group(Cds-Cds(Cds-Cds)H) + group(Cds-CdsOsH) + group(Cdd-CdsCds) + radical(C=COJ) + radical(Cds_S) + radical(C=CJC=C)"""),
)

species(
    label = 'C[C][CH]C#CC=O(25000)',
    structure = SMILES('C[C][CH]C#CC=O'),
    E0 = (590.266,'kJ/mol'),
    modes = [
        HarmonicOscillator(frequencies=([2750,2800,2850,1350,1500,750,1050,1375,1000,2782.5,750,1395,475,1775,1000,2100,2250,500,550,3025,407.5,1350,352.5,291.471,291.606,291.655],'cm^-1')),
        HinderedRotor(inertia=(1.26025,'amu*angstrom^2'), symmetry=1, barrier=(75.9852,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(1.25889,'amu*angstrom^2'), symmetry=1, barrier=(75.9801,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.00198095,'amu*angstrom^2'), symmetry=1, barrier=(0.119627,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(1.25982,'amu*angstrom^2'), symmetry=1, barrier=(75.9808,'kJ/mol'), semiclassical=False),
    ],
    spinMultiplicity = 4,
    opticalIsomers = 1,
    molecularWeight = (93.1033,'amu'),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[1.27032,0.0598299,-6.34887e-05,3.78954e-08,-9.20111e-12,71091.1,24.5361], Tmin=(100,'K'), Tmax=(997.527,'K')), NASAPolynomial(coeffs=[10.2747,0.0237232,-9.19431e-06,1.60941e-09,-1.07119e-13,69294.7,-18.8823], Tmin=(997.527,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(590.266,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(266.063,'J/(mol*K)'), comment="""Thermo group additivity estimation: group(Cs-CsCsHH) + group(Cs-CtCsHH) + group(Cs-CsHHH) + group(Cds-O2d(Cds-Cds)H) + group(Ct-CtCs) + group(Ct-CtCs) + radical(Sec_Propargyl) + radical(CCJ2_triplet)"""),
)

species(
    label = 'C[C]=C[CH][C]=C=O(25001)',
    structure = SMILES('C[C]=C[CH][C]=C=O'),
    E0 = (507.494,'kJ/mol'),
    modes = [
        HarmonicOscillator(frequencies=([1670,1700,300,440,2750,2800,2850,1350,1500,750,1050,1375,1000,2120,512.5,787.5,3010,987.5,1337.5,450,1655,3025,407.5,1350,352.5,252.04,252.718],'cm^-1')),
        HinderedRotor(inertia=(0.607698,'amu*angstrom^2'), symmetry=1, barrier=(27.8129,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.387866,'amu*angstrom^2'), symmetry=1, barrier=(17.4092,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(1.56716,'amu*angstrom^2'), symmetry=1, barrier=(71.7093,'kJ/mol'), semiclassical=False),
    ],
    spinMultiplicity = 4,
    opticalIsomers = 1,
    molecularWeight = (93.1033,'amu'),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[1.39749,0.0588972,-5.86907e-05,3.22675e-08,-7.35398e-12,61129.8,21.9414], Tmin=(100,'K'), Tmax=(1044.08,'K')), NASAPolynomial(coeffs=[9.95328,0.0261183,-1.15974e-05,2.19699e-09,-1.536e-13,59343.2,-19.7042], Tmin=(1044.08,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(507.494,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(270.22,'J/(mol*K)'), comment="""Thermo group additivity estimation: group(Cs-Cd(CCO)HH) + group(Cs-(Cds-Cds)HHH) + group(Cds-CdsCsH) + group(Cds-CdsCsH) + group(Cds-(Cdd-O2d)CsH) + radical(C=CCJC=C=O) + radical(CCCJ=C=O) + radical(Cds_S)"""),
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
    label = '[CH2][C]=CC=C1[CH]O1(25002)',
    structure = SMILES('[CH2][C]=CC=C1[CH]O1'),
    E0 = (513.224,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
    molecularWeight = (93.1033,'amu'),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[1.41469,0.0316685,5.96199e-05,-1.18812e-07,5.33034e-11,61843.1,21.8233], Tmin=(100,'K'), Tmax=(918.338,'K')), NASAPolynomial(coeffs=[23.7004,0.000106955,4.17146e-06,-8.82532e-10,5.31422e-14,54987.6,-98.8338], Tmin=(918.338,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(513.224,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(274.378,'J/(mol*K)'), comment="""Thermo group additivity estimation: group(O2s-Cs(Cds-Cd)) + group(Cs-(Cds-Cds)OsHH) + group(Cs-(Cds-Cds)HHH) + group(Cds-CdsCsOs) + group(Cds-CdsCsH) + group(Cds-Cds(Cds-Cds)H) + group(Cds-Cds(Cds-Cds)H) + ring(methyleneoxirane) + radical(C=CC=CCJ) + radical(Cds_S) + radical(C=CCJO)"""),
)

species(
    label = '[CH2]C1=C[CH]C1=C[O](25003)',
    structure = SMILES('[CH2]C1=C[CH]C1=C[O]'),
    E0 = (368.028,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
    molecularWeight = (93.1033,'amu'),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[1.50881,0.0231478,9.69524e-05,-1.65488e-07,7.19242e-11,44382.7,19.2772], Tmin=(100,'K'), Tmax=(914.154,'K')), NASAPolynomial(coeffs=[26.7354,-0.0047662,7.43622e-06,-1.52245e-09,9.53615e-14,36324.7,-119.01], Tmin=(914.154,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(368.028,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(278.535,'J/(mol*K)'), comment="""Thermo group additivity estimation: group(O2s-(Cds-Cd)H) + group(Cs-(Cds-Cds)(Cds-Cds)HH) + group(Cs-(Cds-Cds)HHH) + group(Cds-Cds(Cds-Cds)Cs) + group(Cds-Cds(Cds-Cds)Cs) + group(Cds-CdsCsH) + group(Cds-CdsOsH) + ring(Cyclobutene) + radical(C=CCJC=C) + radical(Allyl_P) + radical(C=COJ)"""),
)

species(
    label = '[O]C=C1[CH]C=[C]C1(25004)',
    structure = SMILES('[O]C=C1[CH]C=[C]C1'),
    E0 = (407.272,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
    molecularWeight = (93.1033,'amu'),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[2.2739,0.0133932,9.4022e-05,-1.39397e-07,5.65188e-11,49068.3,21.6114], Tmin=(100,'K'), Tmax=(934.705,'K')), NASAPolynomial(coeffs=[17.8005,0.0097301,-8.50691e-07,1.28969e-10,-1.92855e-14,43423.2,-66.9176], Tmin=(934.705,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(407.272,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(282.692,'J/(mol*K)'), comment="""Thermo group additivity estimation: group(O2s-(Cds-Cd)H) + group(Cs-(Cds-Cds)(Cds-Cds)HH) + group(Cs-(Cds-Cds)(Cds-Cds)HH) + group(Cds-CdsCsCs) + group(Cds-CdsCsH) + group(Cds-CdsCsH) + group(Cds-CdsOsH) + ring(4-Methylenecyclopentene) + radical(C=CCJC=C) + radical(C=COJ) + radical(cyclopentene-vinyl)"""),
)

species(
    label = '[CH2]C1=CC=[C]C1[O](24930)',
    structure = SMILES('[CH2]C1=CC=[C]C1[O]'),
    E0 = (524.654,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
    molecularWeight = (93.1033,'amu'),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[1.873,0.0307504,3.43052e-05,-6.91874e-08,2.91809e-11,63192.3,21.7621], Tmin=(100,'K'), Tmax=(964.422,'K')), NASAPolynomial(coeffs=[14.7831,0.0168419,-5.71063e-06,1.08901e-09,-8.30957e-14,58858.9,-49.6105], Tmin=(964.422,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(524.654,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(278.535,'J/(mol*K)'), comment="""Thermo group additivity estimation: group(O2s-CsH) + group(Cs-(Cds-Cds)(Cds-Cds)OsH) + group(Cs-(Cds-Cds)HHH) + group(Cds-CdsCsCs) + group(Cds-CdsCsH) + group(Cds-Cds(Cds-Cds)H) + group(Cds-Cds(Cds-Cds)H) + ring(Cyclopentadiene) + radical(CC(C)OJ) + radical(1,3-cyclopentadiene-vinyl-1) + radical(C=CC=CCJ)"""),
)

species(
    label = '[O]C1[C]=CC=[C]C1(24844)',
    structure = SMILES('[O]C1[C]=CC=[C]C1'),
    E0 = (621.61,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
    molecularWeight = (93.1033,'amu'),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[1.88485,0.0328097,2.13486e-05,-5.20611e-08,2.21927e-11,74850.8,20.885], Tmin=(100,'K'), Tmax=(983.273,'K')), NASAPolynomial(coeffs=[13.6533,0.0180947,-6.78906e-06,1.31394e-09,-9.83113e-14,70933.5,-43.8438], Tmin=(983.273,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(621.61,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(282.692,'J/(mol*K)'), comment="""Thermo group additivity estimation: group(O2s-CsH) + group(Cs-(Cds-Cds)CsOsH) + group(Cs-(Cds-Cds)CsHH) + group(Cds-CdsCsH) + group(Cds-CdsCsH) + group(Cds-Cds(Cds-Cds)H) + group(Cds-Cds(Cds-Cds)H) + ring(1,3-Cyclohexadiene) + radical(CC(C)OJ) + radical(Cds_S) + radical(Cds_S)"""),
)

species(
    label = '[CH2]C=CC=C=C=O(25005)',
    structure = SMILES('[CH2]C=CC=C=C=O'),
    E0 = (293.297,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
    molecularWeight = (93.1033,'amu'),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[2.33764,0.0243379,2.66865e-05,-5.49926e-08,2.3737e-11,35346.4,6.65685], Tmin=(100,'K'), Tmax=(945.854,'K')), NASAPolynomial(coeffs=[12.4055,0.0125298,-3.38231e-06,5.93006e-10,-4.52178e-14,32065.6,-48.6299], Tmin=(945.854,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(293.297,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(274.378,'J/(mol*K)'), comment="""Thermo group additivity estimation: group(Cs-(Cds-Cds)HHH) + group(Cds-CdsCsH) + group(Cds-Cds(Cds-Cds)H) + group(Cds-Cds(Cds-Cds)H) + group(Cdd-CdsCds) + radical(C=CC=CCJ)"""),
)

species(
    label = '[CH2][C]=[C]C[C]=C[O](25006)',
    structure = SMILES('[CH2][C]=[C]C[C]=C[O]'),
    E0 = (850.419,'kJ/mol'),
    modes = [
        HarmonicOscillator(frequencies=([1670,1685,1700,300,370,440,2750,2850,1437.5,1250,1305,750,350,3010,987.5,1337.5,450,1655,3000,3100,440,815,1455,1000,347.806,348.061,348.464],'cm^-1')),
        HinderedRotor(inertia=(0.223581,'amu*angstrom^2'), symmetry=1, barrier=(19.1664,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.223594,'amu*angstrom^2'), symmetry=1, barrier=(19.1769,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.222263,'amu*angstrom^2'), symmetry=1, barrier=(19.1656,'kJ/mol'), semiclassical=False),
    ],
    spinMultiplicity = 6,
    opticalIsomers = 1,
    molecularWeight = (93.1033,'amu'),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[0.804933,0.0644284,-6.5943e-05,3.40461e-08,-6.87889e-12,102402,29.0695], Tmin=(100,'K'), Tmax=(1211.68,'K')), NASAPolynomial(coeffs=[15.4795,0.0159845,-5.97184e-06,1.04993e-09,-7.09389e-14,98845.6,-44.5445], Tmin=(1211.68,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(850.419,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(270.22,'J/(mol*K)'), comment="""Thermo group additivity estimation: group(O2s-(Cds-Cd)H) + group(Cs-(Cds-Cds)(Cds-Cds)HH) + group(Cs-(Cds-Cds)HHH) + group(Cds-CdsCsH) + group(Cds-CdsCsH) + group(Cds-CdsCsH) + group(Cds-CdsOsH) + radical(C=COJ) + radical(Allyl_P) + radical(Cds_S) + radical(Cds_S) + radical(Cds_S)"""),
)

species(
    label = '[CH2][C]=C[C]=C[CH][O](25007)',
    structure = SMILES('[CH2][C]=C[C]=C[CH][O]'),
    E0 = (766.628,'kJ/mol'),
    modes = [
        HarmonicOscillator(frequencies=([1670,1700,300,440,3025,407.5,1350,352.5,2995,3025,975,1000,1300,1375,400,500,1630,1680,3000,3100,440,815,1455,1000,440.523,440.771,441.197],'cm^-1')),
        HinderedRotor(inertia=(0.000864431,'amu*angstrom^2'), symmetry=1, barrier=(0.119627,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.385383,'amu*angstrom^2'), symmetry=1, barrier=(53.0998,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.38448,'amu*angstrom^2'), symmetry=1, barrier=(53.0784,'kJ/mol'), semiclassical=False),
    ],
    spinMultiplicity = 6,
    opticalIsomers = 1,
    molecularWeight = (93.1033,'amu'),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[1.46312,0.052414,-4.18821e-05,1.82756e-08,-3.30526e-12,92298.1,27.5545], Tmin=(100,'K'), Tmax=(1298.77,'K')), NASAPolynomial(coeffs=[10.5444,0.024445,-9.57972e-06,1.69461e-09,-1.13579e-13,89939.2,-18.6315], Tmin=(1298.77,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(766.628,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(270.22,'J/(mol*K)'), comment="""Thermo group additivity estimation: group(O2s-CsH) + group(Cs-(Cds-Cds)OsHH) + group(Cs-(Cds-Cds)HHH) + group(Cds-CdsCsH) + group(Cds-CdsCsH) + group(Cds-Cds(Cds-Cds)H) + group(Cds-Cds(Cds-Cds)H) + radical(C=CCJO) + radical(C=CC=CCJ) + radical(CCOJ) + radical(Cds_S) + radical(C=CJC=C)"""),
)

species(
    label = '[CH2][C]=C[CH]C=[C][O](25008)',
    structure = SMILES('[CH2][C]=C[CH]C=[C][O]'),
    E0 = (716.204,'kJ/mol'),
    modes = [
        HarmonicOscillator(frequencies=([1670,1700,300,440,3025,407.5,1350,352.5,2995,3025,975,1000,1300,1375,400,500,1630,1680,3000,3100,440,815,1455,1000,600.337,600.455,601.309],'cm^-1')),
        HinderedRotor(inertia=(0.0862025,'amu*angstrom^2'), symmetry=1, barrier=(22.0797,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.348912,'amu*angstrom^2'), symmetry=1, barrier=(89.3139,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.349232,'amu*angstrom^2'), symmetry=1, barrier=(89.3124,'kJ/mol'), semiclassical=False),
    ],
    spinMultiplicity = 6,
    opticalIsomers = 1,
    molecularWeight = (93.1033,'amu'),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[1.45992,0.0439142,-3.40976e-06,-3.29505e-08,1.75507e-11,86241.8,28.569], Tmin=(100,'K'), Tmax=(946.739,'K')), NASAPolynomial(coeffs=[15.2739,0.0150819,-4.51828e-06,7.78354e-10,-5.63791e-14,82302.7,-44.3091], Tmin=(946.739,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(716.204,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(270.22,'J/(mol*K)'), comment="""Thermo group additivity estimation: group(O2s-(Cds-Cd)H) + group(Cs-(Cds-Cds)(Cds-Cds)HH) + group(Cs-(Cds-Cds)HHH) + group(Cds-CdsCsH) + group(Cds-CdsCsH) + group(Cds-CdsCsH) + group(Cds-CdsOsH) + radical(Cds_S) + radical(C=CCJC=C) + radical(Allyl_P) + radical(C=COJ) + radical(C=CJO)"""),
)

species(
    label = '[CH2][C]=[C][CH]C=C[O](25009)',
    structure = SMILES('[CH2][C]=[C][CH]C=C[O]'),
    E0 = (714.302,'kJ/mol'),
    modes = [
        HarmonicOscillator(frequencies=([1670,1700,300,440,3025,407.5,1350,352.5,2995,3025,975,1000,1300,1375,400,500,1630,1680,3000,3100,440,815,1455,1000,463.822,464.06,464.08],'cm^-1')),
        HinderedRotor(inertia=(0.21979,'amu*angstrom^2'), symmetry=1, barrier=(33.5841,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.219961,'amu*angstrom^2'), symmetry=1, barrier=(33.5866,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.219802,'amu*angstrom^2'), symmetry=1, barrier=(33.5881,'kJ/mol'), semiclassical=False),
    ],
    spinMultiplicity = 6,
    opticalIsomers = 1,
    molecularWeight = (93.1033,'amu'),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[1.22352,0.0466796,-1.04833e-06,-4.28059e-08,2.29479e-11,86023.9,26.7627], Tmin=(100,'K'), Tmax=(930.61,'K')), NASAPolynomial(coeffs=[18.1793,0.0105324,-1.99271e-06,2.85902e-10,-2.28299e-14,81277.5,-62.3658], Tmin=(930.61,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(714.302,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(270.22,'J/(mol*K)'), comment="""Thermo group additivity estimation: group(O2s-(Cds-Cd)H) + group(Cs-(Cds-Cds)(Cds-Cds)HH) + group(Cs-(Cds-Cds)HHH) + group(Cds-CdsCsH) + group(Cds-CdsCsH) + group(Cds-CdsCsH) + group(Cds-CdsOsH) + radical(Cds_S) + radical(Cds_S) + radical(C=COJ) + radical(C=CCJC=C) + radical(Allyl_P)"""),
)

species(
    label = '[CH2][C]=[C]C=[C]C[O](25010)',
    structure = SMILES('[CH2][C]=[C]C=[C]C[O]'),
    E0 = (887.174,'kJ/mol'),
    modes = [
        HarmonicOscillator(frequencies=([1670,1685,1700,300,370,440,2750,2850,1437.5,1250,1305,750,350,3010,987.5,1337.5,450,1655,3000,3100,440,815,1455,1000,180,180,690.912],'cm^-1')),
        HinderedRotor(inertia=(0.0821643,'amu*angstrom^2'), symmetry=1, barrier=(1.88912,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.0956578,'amu*angstrom^2'), symmetry=1, barrier=(2.19936,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(2.5116,'amu*angstrom^2'), symmetry=1, barrier=(57.7467,'kJ/mol'), semiclassical=False),
    ],
    spinMultiplicity = 6,
    opticalIsomers = 1,
    molecularWeight = (93.1033,'amu'),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[1.33467,0.0667767,-0.000102424,9.33456e-08,-3.29722e-11,106791,27.8915], Tmin=(100,'K'), Tmax=(876.197,'K')), NASAPolynomial(coeffs=[4.2578,0.0347154,-1.54945e-05,2.82486e-09,-1.88392e-13,106997,18.2753], Tmin=(876.197,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(887.174,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(270.22,'J/(mol*K)'), comment="""Thermo group additivity estimation: group(O2s-CsH) + group(Cs-(Cds-Cds)OsHH) + group(Cs-(Cds-Cds)HHH) + group(Cds-CdsCsH) + group(Cds-CdsCsH) + group(Cds-Cds(Cds-Cds)H) + group(Cds-Cds(Cds-Cds)H) + radical(Cds_S) + radical(CCOJ) + radical(C=CC=CCJ) + radical(C=CJC=C) + radical(Cds_S)"""),
)

species(
    label = 'C=[C]C[C]=C=C[O](24306)',
    structure = SMILES('C=[C]C[C]=C=C[O]'),
    E0 = (637.682,'kJ/mol'),
    modes = [
        HarmonicOscillator(frequencies=([1670,1700,300,440,2750,2850,1437.5,1250,1305,750,350,540,610,2055,3010,987.5,1337.5,450,1655,2950,3100,1380,975,1025,1650,180,180,180],'cm^-1')),
        HinderedRotor(inertia=(0.996479,'amu*angstrom^2'), symmetry=1, barrier=(22.911,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.99986,'amu*angstrom^2'), symmetry=1, barrier=(22.9888,'kJ/mol'), semiclassical=False),
    ],
    spinMultiplicity = 4,
    opticalIsomers = 1,
    molecularWeight = (93.1033,'amu'),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[0.788263,0.0629847,-6.2157e-05,3.07557e-08,-5.93811e-12,76817.5,27.4308], Tmin=(100,'K'), Tmax=(1268.94,'K')), NASAPolynomial(coeffs=[16.0234,0.0149606,-5.38889e-06,9.31594e-10,-6.23959e-14,72951,-49.6988], Tmin=(1268.94,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(637.682,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(274.378,'J/(mol*K)'), comment="""Thermo group additivity estimation: group(O2s-(Cds-Cd)H) + group(Cs-(Cds-Cds)(Cds-Cds)HH) + group(Cds-CdsCsH) + group(Cds-CdsCsH) + group(Cds-CdsHH) + group(Cds-CdsOsH) + group(Cdd-CdsCds) + radical(C=COJ) + radical(Cds_S) + radical(Cds_S)"""),
)

species(
    label = '[CH]=[C]CC=C=C[O](22626)',
    structure = SMILES('[CH]=[C]CC=C=C[O]'),
    E0 = (646.936,'kJ/mol'),
    modes = [
        HarmonicOscillator(frequencies=([1685,370,2750,2850,1437.5,1250,1305,750,350,2995,3025,975,1000,1300,1375,400,500,1630,1680,3120,650,792.5,1650,540,610,2055,180,180],'cm^-1')),
        HinderedRotor(inertia=(1.06585,'amu*angstrom^2'), symmetry=1, barrier=(24.506,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(1.06786,'amu*angstrom^2'), symmetry=1, barrier=(24.5523,'kJ/mol'), semiclassical=False),
    ],
    spinMultiplicity = 4,
    opticalIsomers = 1,
    molecularWeight = (93.1033,'amu'),
    collisionModel = TransportData(shapeIndex=2, epsilon=(3942.94,'J/mol'), sigma=(6.38069,'angstroms'), dipoleMoment=(0,'C*m'), polarizability=(0,'angstroms^3'), rotrelaxcollnum=0, comment="""Epsilon & sigma estimated with Tc=615.88 K, Pc=34.44 bar (from Joback method)"""),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[0.874198,0.0604722,-4.9312e-05,1.32472e-08,1.2772e-12,77928.4,27.0485], Tmin=(100,'K'), Tmax=(985.722,'K')), NASAPolynomial(coeffs=[16.5106,0.0141667,-4.93912e-06,8.83067e-10,-6.24589e-14,74012.8,-52.3884], Tmin=(985.722,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(646.936,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(274.378,'J/(mol*K)'), comment="""Thermo group additivity estimation: group(O2s-(Cds-Cd)H) + group(Cs-(Cds-Cds)(Cds-Cds)HH) + group(Cds-CdsCsH) + group(Cds-CdsCsH) + group(Cds-CdsHH) + group(Cds-CdsOsH) + group(Cdd-CdsCds) + radical(Cds_S) + radical(C=COJ) + radical(Cds_P)"""),
)

species(
    label = 'C=[C]C[CH][C]=C=O(24307)',
    structure = SMILES('C=[C]C[CH][C]=C=O'),
    E0 = (565.963,'kJ/mol'),
    modes = [
        HarmonicOscillator(frequencies=([1670,1700,300,440,2750,2850,1437.5,1250,1305,750,350,2950,3100,1380,975,1025,1650,2120,512.5,787.5,3025,407.5,1350,352.5,304.114,305.404,1203.3],'cm^-1')),
        HinderedRotor(inertia=(0.159997,'amu*angstrom^2'), symmetry=1, barrier=(10.5335,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.886339,'amu*angstrom^2'), symmetry=1, barrier=(57.9697,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.307049,'amu*angstrom^2'), symmetry=1, barrier=(20.1994,'kJ/mol'), semiclassical=False),
    ],
    spinMultiplicity = 4,
    opticalIsomers = 1,
    molecularWeight = (93.1033,'amu'),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[1.42967,0.0609813,-7.03206e-05,4.80318e-08,-1.39051e-11,68158.3,25.9364], Tmin=(100,'K'), Tmax=(825.376,'K')), NASAPolynomial(coeffs=[7.79175,0.0301483,-1.42851e-05,2.7704e-09,-1.95504e-13,67108.1,-3.53579], Tmin=(825.376,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(565.963,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(270.22,'J/(mol*K)'), comment="""Thermo group additivity estimation: group(Cs-(Cds-Cds)CsHH) + group(Cs-(Cds-Cdd-O2d)CsHH) + group(Cds-CdsCsH) + group(Cds-(Cdd-O2d)CsH) + group(Cds-CdsHH) + radical(Cds_S) + radical(CCJC(C)=C=O) + radical(CCCJ=C=O)"""),
)

species(
    label = 'C=[C]C1C=[C]C1[O](25011)',
    structure = SMILES('C=[C]C1C=[C]C1[O]'),
    E0 = (764.051,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
    molecularWeight = (93.1033,'amu'),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[1.41886,0.0487991,-2.9315e-05,2.85123e-09,2.30414e-12,91993.8,23.5359], Tmin=(100,'K'), Tmax=(1081.73,'K')), NASAPolynomial(coeffs=[13.0213,0.019673,-8.03115e-06,1.50808e-09,-1.06548e-13,88677.6,-37.0763], Tmin=(1081.73,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(764.051,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(278.535,'J/(mol*K)'), comment="""Thermo group additivity estimation: group(O2s-CsH) + group(Cs-(Cds-Cds)(Cds-Cds)CsH) + group(Cs-(Cds-Cds)CsOsH) + group(Cds-CdsCsH) + group(Cds-CdsCsH) + group(Cds-CdsCsH) + group(Cds-CdsHH) + ring(Cyclobutene) + radical(CC(C)OJ) + radical(Cds_S) + radical(cyclobutene-vinyl)"""),
)

species(
    label = 'C=[C]CC=C=C=O(24300)',
    structure = SMILES('C=[C]CC=C=C=O'),
    E0 = (443.937,'kJ/mol'),
    modes = [
        HarmonicOscillator(frequencies=([2950,3100,1380,975,1025,1650,540,610,2055,2750,2850,1437.5,1250,1305,750,350,2120,512.5,787.5,1685,370,3010,987.5,1337.5,450,1655,4000,4000],'cm^-1')),
        HinderedRotor(inertia=(0.00088996,'amu*angstrom^2'), symmetry=1, barrier=(10.1047,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.000888337,'amu*angstrom^2'), symmetry=1, barrier=(10.0862,'kJ/mol'), semiclassical=False),
    ],
    spinMultiplicity = 2,
    opticalIsomers = 1,
    molecularWeight = (93.1033,'amu'),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[2.2116,0.0334925,-1.42842e-05,-2.85593e-09,2.74141e-12,53462.5,9.05176], Tmin=(100,'K'), Tmax=(1115.95,'K')), NASAPolynomial(coeffs=[9.57185,0.0176473,-7.14894e-06,1.34229e-09,-9.45158e-14,51163.7,-30.2043], Tmin=(1115.95,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(443.937,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(274.378,'J/(mol*K)'), comment="""Thermo group additivity estimation: group(Cs-(Cds-Cds)(Cds-Cds)HH) + group(Cds-CdsCsH) + group(Cds-CdsCsH) + group(Cds-CdsHH) + group(Cdd-CdsCds) + radical(Cds_S)"""),
)

species(
    label = '[CH2]C#CC[C]=C[O](25012)',
    structure = SMILES('[CH2]C#CC[C]=C[O]'),
    E0 = (542.868,'kJ/mol'),
    modes = [
        HarmonicOscillator(frequencies=([1685,370,2750,2850,1437.5,1250,1305,750,350,3010,987.5,1337.5,450,1655,2100,2250,500,550,3000,3100,440,815,1455,1000,260.722,260.724,260.725],'cm^-1')),
        HinderedRotor(inertia=(0.638997,'amu*angstrom^2'), symmetry=1, barrier=(30.8235,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.00248,'amu*angstrom^2'), symmetry=1, barrier=(0.119627,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.638997,'amu*angstrom^2'), symmetry=1, barrier=(30.8235,'kJ/mol'), semiclassical=False),
    ],
    spinMultiplicity = 4,
    opticalIsomers = 1,
    molecularWeight = (93.1033,'amu'),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[1.22495,0.0523846,-3.1887e-05,-1.06973e-09,5.55251e-12,65399.5,26.1705], Tmin=(100,'K'), Tmax=(978.736,'K')), NASAPolynomial(coeffs=[14.752,0.0161948,-5.68602e-06,1.01607e-09,-7.16801e-14,61837.1,-43.4708], Tmin=(978.736,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(542.868,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(270.22,'J/(mol*K)'), comment="""Thermo group additivity estimation: group(O2s-(Cds-Cd)H) + group(Cs-(Cds-Cds)CtHH) + group(Cs-CtHHH) + group(Cds-CdsCsH) + group(Cds-CdsOsH) + group(Ct-CtCs) + group(Ct-CtCs) + radical(Propargyl) + radical(Cds_S) + radical(C=COJ)"""),
)

species(
    label = 'C=C=CC[C][C]=O(25013)',
    structure = SMILES('C=C=CC[C][C]=O'),
    E0 = (625.653,'kJ/mol'),
    modes = [
        HarmonicOscillator(frequencies=([1855,455,950,2750,2850,1437.5,1250,1305,750,350,540,610,2055,3010,987.5,1337.5,450,1655,2950,3100,1380,975,1025,1650,180.608,180.65,180.67],'cm^-1')),
        HinderedRotor(inertia=(0.542468,'amu*angstrom^2'), symmetry=1, barrier=(12.5464,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.541675,'amu*angstrom^2'), symmetry=1, barrier=(12.5462,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.966974,'amu*angstrom^2'), symmetry=1, barrier=(22.3962,'kJ/mol'), semiclassical=False),
    ],
    spinMultiplicity = 4,
    opticalIsomers = 1,
    molecularWeight = (93.1033,'amu'),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[1.20057,0.0655387,-8.22355e-05,5.84425e-08,-1.70471e-11,75345.9,25.997], Tmin=(100,'K'), Tmax=(830.499,'K')), NASAPolynomial(coeffs=[9.23153,0.0268579,-1.23715e-05,2.35968e-09,-1.64596e-13,74012,-11.256], Tmin=(830.499,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(625.653,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(270.22,'J/(mol*K)'), comment="""Thermo group additivity estimation: group(Cs-(Cds-Cds)CsHH) + group(Cs-(Cds-O2d)CsHH) + group(Cds-CdsCsH) + group(Cds-OdCsH) + group(Cds-CdsHH) + group(Cdd-CdsCds) + radical(CCCJ=O) + radical(CCJ2_triplet)"""),
)

species(
    label = 'C#C[CH]C[C]=C[O](25014)',
    structure = SMILES('C#C[CH]C[C]=C[O]'),
    E0 = (546.525,'kJ/mol'),
    modes = [
        HarmonicOscillator(frequencies=([750,770,3400,2100,3025,407.5,1350,352.5,2175,525,2750,2850,1437.5,1250,1305,750,350,1685,370,3010,987.5,1337.5,450,1655,296.018,296.027,296.039],'cm^-1')),
        HinderedRotor(inertia=(1.15081,'amu*angstrom^2'), symmetry=1, barrier=(71.5537,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.370722,'amu*angstrom^2'), symmetry=1, barrier=(23.0554,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(1.15049,'amu*angstrom^2'), symmetry=1, barrier=(71.5538,'kJ/mol'), semiclassical=False),
    ],
    spinMultiplicity = 4,
    opticalIsomers = 1,
    molecularWeight = (93.1033,'amu'),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[0.401384,0.0667083,-7.02614e-05,3.71121e-08,-7.3853e-12,65872,27.3776], Tmin=(100,'K'), Tmax=(1409.93,'K')), NASAPolynomial(coeffs=[16.6854,0.0113457,-1.61232e-06,4.21859e-11,5.14822e-15,62191.1,-53.5474], Tmin=(1409.93,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(546.525,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(270.22,'J/(mol*K)'), comment="""Thermo group additivity estimation: group(O2s-(Cds-Cd)H) + group(Cs-(Cds-Cds)CsHH) + group(Cs-CtCsHH) + group(Cds-CdsCsH) + group(Cds-CdsOsH) + group(Ct-CtCs) + group(Ct-CtH) + radical(Cds_S) + radical(Sec_Propargyl) + radical(C=COJ)"""),
)

species(
    label = 'C=C1[CH]C=[C][CH]O1(24914)',
    structure = SMILES('C=C1[CH]C=[C][CH]O1'),
    E0 = (389.249,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
    molecularWeight = (93.1033,'amu'),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[2.00434,0.021333,7.53505e-05,-1.12198e-07,4.31238e-11,46907.8,16.018], Tmin=(100,'K'), Tmax=(991.103,'K')), NASAPolynomial(coeffs=[16.225,0.0204127,-8.72692e-06,1.84816e-09,-1.45646e-13,41315.3,-66.4543], Tmin=(991.103,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(389.249,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(282.692,'J/(mol*K)'), comment="""Thermo group additivity estimation: group(O2s-Cs(Cds-Cd)) + group(Cs-(Cds-Cds)(Cds-Cds)HH) + group(Cs-(Cds-Cds)OsHH) + group(Cds-CdsCsOs) + group(Cds-CdsCsH) + group(Cds-CdsCsH) + group(Cds-CdsHH) + ring(Cyclohexane) + radical(C=CCJ(O)C) + radical(Cds_S) + radical(C=CCJCO)"""),
)

species(
    label = '[C]1=C[CH][C]=COC1(24829)',
    structure = SMILES('[C]1=C[CH][C]=COC1'),
    E0 = (536.895,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
    molecularWeight = (93.1033,'amu'),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[2.15374,-0.0253958,0.00028947,-4.13616e-07,1.74365e-10,64702.4,21.0251], Tmin=(100,'K'), Tmax=(895.873,'K')), NASAPolynomial(coeffs=[45.6272,-0.0459226,3.32084e-05,-6.64532e-09,4.44299e-13,49947.4,-222.805], Tmin=(895.873,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(536.895,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(282.692,'J/(mol*K)'), comment="""Thermo group additivity estimation: group(O2s-Cs(Cds-Cd)) + group(Cs-(Cds-Cds)(Cds-Cds)HH) + group(Cs-(Cds-Cds)OsHH) + group(Cds-CdsCsH) + group(Cds-CdsCsH) + group(Cds-CdsCsH) + group(Cds-CdsOsH) + ring(Cycloheptane) + radical(C=CCJC=C) + radical(Cds_S) + radical(Cds_S)"""),
)

species(
    label = '[CH2]C=[C][CH][C]=C[O](25015)',
    structure = SMILES('[CH2]C=[C][CH][C]=C[O]'),
    E0 = (714.302,'kJ/mol'),
    modes = [
        HarmonicOscillator(frequencies=([1670,1700,300,440,3025,407.5,1350,352.5,2995,3025,975,1000,1300,1375,400,500,1630,1680,3000,3100,440,815,1455,1000,463.944,464.177,464.269],'cm^-1')),
        HinderedRotor(inertia=(0.219077,'amu*angstrom^2'), symmetry=1, barrier=(33.5779,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.219855,'amu*angstrom^2'), symmetry=1, barrier=(33.5739,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.221035,'amu*angstrom^2'), symmetry=1, barrier=(33.6071,'kJ/mol'), semiclassical=False),
    ],
    spinMultiplicity = 6,
    opticalIsomers = 1,
    molecularWeight = (93.1033,'amu'),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[1.22354,0.0466793,-1.04735e-06,-4.28072e-08,2.29485e-11,86023.9,26.7627], Tmin=(100,'K'), Tmax=(930.607,'K')), NASAPolynomial(coeffs=[18.1792,0.0105325,-1.99278e-06,2.85918e-10,-2.28312e-14,81277.5,-62.3654], Tmin=(930.607,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(714.302,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(270.22,'J/(mol*K)'), comment="""Thermo group additivity estimation: group(O2s-(Cds-Cd)H) + group(Cs-(Cds-Cds)(Cds-Cds)HH) + group(Cs-(Cds-Cds)HHH) + group(Cds-CdsCsH) + group(Cds-CdsCsH) + group(Cds-CdsCsH) + group(Cds-CdsOsH) + radical(Allyl_P) + radical(Cds_S) + radical(C=CCJC=C) + radical(Cds_S) + radical(C=COJ)"""),
)

species(
    label = '[CH]C=C[CH][C]=C[O](25016)',
    structure = SMILES('[CH]C=C[CH][C]=C[O]'),
    E0 = (695.645,'kJ/mol'),
    modes = [
        HarmonicOscillator(frequencies=([3025,407.5,1350,352.5,2995,3010,3025,975,987.5,1000,1300,1337.5,1375,400,450,500,1630,1655,1680,1685,370,343.938,343.938,343.939,343.939,343.939,343.94],'cm^-1')),
        HinderedRotor(inertia=(0.594291,'amu*angstrom^2'), symmetry=1, barrier=(49.8872,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.59429,'amu*angstrom^2'), symmetry=1, barrier=(49.8872,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.594289,'amu*angstrom^2'), symmetry=1, barrier=(49.8872,'kJ/mol'), semiclassical=False),
    ],
    spinMultiplicity = 6,
    opticalIsomers = 1,
    molecularWeight = (93.1033,'amu'),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[1.25669,0.0444984,1.55777e-05,-5.93077e-08,2.81736e-11,83780.1,26.9982], Tmin=(100,'K'), Tmax=(936.133,'K')), NASAPolynomial(coeffs=[16.9982,0.0172119,-4.75437e-06,7.88044e-10,-5.74018e-14,79081.3,-57.2624], Tmin=(936.133,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(695.645,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(270.22,'J/(mol*K)'), comment="""Thermo group additivity estimation: group(O2s-(Cds-Cd)H) + group(Cs-(Cds-Cds)(Cds-Cds)HH) + group(Cs-(Cds-Cds)HHH) + group(Cds-CdsCsH) + group(Cds-CdsCsH) + group(Cds-CdsCsH) + group(Cds-CdsOsH) + radical(Cds_S) + radical(AllylJ2_triplet) + radical(C=COJ) + radical(C=CCJC=C)"""),
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
    label = '[C]C=C[C]=C(19282)',
    structure = SMILES('[C]C=C[C]=C'),
    E0 = (926.369,'kJ/mol'),
    modes = [
        HarmonicOscillator(frequencies=([2950,3100,1380,975,1025,1650,2995,3025,975,1000,1300,1375,400,500,1630,1680,1685,370,180,180],'cm^-1')),
        HinderedRotor(inertia=(1.54615,'amu*angstrom^2'), symmetry=1, barrier=(35.5489,'kJ/mol'), semiclassical=False),
    ],
    spinMultiplicity = 5,
    opticalIsomers = 1,
    molecularWeight = (64.0853,'amu'),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[1.60723,0.0446451,-4.44308e-05,2.23547e-08,-4.30798e-12,111509,16.4407], Tmin=(100,'K'), Tmax=(1401.68,'K')), NASAPolynomial(coeffs=[12.6996,0.00908062,-2.18752e-06,2.72849e-10,-1.45836e-14,108784,-39.4491], Tmin=(1401.68,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(926.369,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(203.705,'J/(mol*K)'), comment="""Thermo group additivity estimation: group(Cs-(Cds-Cds)HHH) + group(Cds-CdsCsH) + group(Cds-Cds(Cds-Cds)H) + group(Cds-Cds(Cds-Cds)H) + group(Cds-CdsHH) + radical(C=CJC=C) + radical(CJ3)"""),
)

species(
    label = 'C=[C]C=C=CC=O(25017)',
    structure = SMILES('C=[C]C=C=CC=O'),
    E0 = (310.007,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
    molecularWeight = (93.1033,'amu'),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[1.36062,0.0613328,-6.32705e-05,3.55083e-08,-8.27983e-12,37377.6,21.013], Tmin=(100,'K'), Tmax=(1017.96,'K')), NASAPolynomial(coeffs=[10.0576,0.0271579,-1.29118e-05,2.52752e-09,-1.79961e-13,35607,-21.0996], Tmin=(1017.96,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(310.007,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(274.378,'J/(mol*K)'), comment="""Thermo group additivity estimation: group(Cds-Cds(Cds-Cds)H) + group(Cds-Cds(Cds-Cds)H) + group(Cd-Cd(CO)H) + group(Cds-CdsHH) + group(Cds-O2d(Cds-Cds)H) + group(Cdd-CdsCds) + radical(C=CJC=C)"""),
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
    label = '[CH2][C]=CC#CC[O](25018)',
    structure = SMILES('[CH2][C]=CC#CC[O]'),
    E0 = (624.328,'kJ/mol'),
    modes = [
        HarmonicOscillator(frequencies=([1685,370,2750,2850,1437.5,1250,1305,750,350,3010,987.5,1337.5,450,1655,2100,2250,500,550,3000,3100,440,815,1455,1000,180,180,1436.1],'cm^-1')),
        HinderedRotor(inertia=(0.190746,'amu*angstrom^2'), symmetry=1, barrier=(4.38562,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.0419712,'amu*angstrom^2'), symmetry=1, barrier=(61.6268,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(2.67204,'amu*angstrom^2'), symmetry=1, barrier=(61.4355,'kJ/mol'), semiclassical=False),
    ],
    spinMultiplicity = 4,
    opticalIsomers = 1,
    molecularWeight = (93.1033,'amu'),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[1.72113,0.0543448,-6.62032e-05,5.44939e-08,-1.87392e-11,75167.3,25.0363], Tmin=(100,'K'), Tmax=(844.154,'K')), NASAPolynomial(coeffs=[4.67349,0.0326843,-1.40837e-05,2.56818e-09,-1.73173e-13,74942.2,12.9119], Tmin=(844.154,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(624.328,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(270.22,'J/(mol*K)'), comment="""Thermo group additivity estimation: group(O2s-CsH) + group(Cs-(Cds-Cds)HHH) + group(Cs-CtOsHH) + group(Cds-CdsCsH) + group(Cds-CdsCtH) + group(Ct-CtCs) + group(Ct-Ct(Cds-Cds)) + radical(Cds_S) + radical(CTCC=CCJ) + radical(CCOJ)"""),
)

species(
    label = '[CH2]C#CC=[C]C[O](25019)',
    structure = SMILES('[CH2]C#CC=[C]C[O]'),
    E0 = (659.319,'kJ/mol'),
    modes = [
        HarmonicOscillator(frequencies=([1685,370,2750,2850,1437.5,1250,1305,750,350,3010,987.5,1337.5,450,1655,2100,2250,500,550,3000,3100,440,815,1455,1000,223.053,223.075,1964.71],'cm^-1')),
        HinderedRotor(inertia=(0.393611,'amu*angstrom^2'), symmetry=1, barrier=(13.8966,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(1.0671,'amu*angstrom^2'), symmetry=1, barrier=(37.6754,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.0137548,'amu*angstrom^2'), symmetry=1, barrier=(37.6752,'kJ/mol'), semiclassical=False),
    ],
    spinMultiplicity = 4,
    opticalIsomers = 1,
    molecularWeight = (93.1033,'amu'),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[1.7336,0.0545847,-6.57966e-05,5.5818e-08,-2.0918e-11,79374.8,26.3982], Tmin=(100,'K'), Tmax=(733.052,'K')), NASAPolynomial(coeffs=[4.49405,0.03514,-1.70418e-05,3.32398e-09,-2.34517e-13,79087.9,14.7408], Tmin=(733.052,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(659.319,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(270.22,'J/(mol*K)'), comment="""Thermo group additivity estimation: group(O2s-CsH) + group(Cs-(Cds-Cds)OsHH) + group(Cs-CtHHH) + group(Cds-CdsCsH) + group(Cds-CdsCtH) + group(Ct-CtCs) + group(Ct-Ct(Cds-Cds)) + radical(Cds_S) + radical(CCOJ) + radical(Propargyl)"""),
)

species(
    label = '[CH]=C=CC=[C]C[O](25020)',
    structure = SMILES('[CH]=C=CC=[C]C[O]'),
    E0 = (663.362,'kJ/mol'),
    modes = [
        HarmonicOscillator(frequencies=([3120,650,792.5,1650,2750,2850,1437.5,1250,1305,750,350,540,610,2055,2995,3025,975,1000,1300,1375,400,500,1630,1680,1685,370,180,180],'cm^-1')),
        HinderedRotor(inertia=(0.535255,'amu*angstrom^2'), symmetry=1, barrier=(12.3066,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.716921,'amu*angstrom^2'), symmetry=1, barrier=(16.4834,'kJ/mol'), semiclassical=False),
    ],
    spinMultiplicity = 4,
    opticalIsomers = 1,
    molecularWeight = (93.1033,'amu'),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[1.2322,0.0649858,-8.79826e-05,7.01803e-08,-2.26544e-11,79879.8,26.3124], Tmin=(100,'K'), Tmax=(845.716,'K')), NASAPolynomial(coeffs=[7.81119,0.0280345,-1.20957e-05,2.20216e-09,-1.481e-13,78975.6,-3.09135], Tmin=(845.716,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(663.362,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(274.378,'J/(mol*K)'), comment="""Thermo group additivity estimation: group(O2s-CsH) + group(Cs-(Cds-Cds)OsHH) + group(Cds-CdsCsH) + group(Cds-Cds(Cds-Cds)H) + group(Cds-Cds(Cds-Cds)H) + group(Cds-CdsHH) + group(Cdd-CdsCds) + radical(Cds_S) + radical(C=C=CJ) + radical(CCOJ)"""),
)

species(
    label = 'C=[C]C[C][C]=C[O](25021)',
    structure = SMILES('C=[C]C[C][C]=C[O]'),
    E0 = (920.269,'kJ/mol'),
    modes = [
        HarmonicOscillator(frequencies=([1670,1700,300,440,2750,2850,1437.5,1250,1305,750,350,3010,987.5,1337.5,450,1655,2950,3100,1380,975,1025,1650,297.965,297.965,297.965,297.965,297.965],'cm^-1')),
        HinderedRotor(inertia=(0.333656,'amu*angstrom^2'), symmetry=1, barrier=(21.0212,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.333656,'amu*angstrom^2'), symmetry=1, barrier=(21.0212,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.333656,'amu*angstrom^2'), symmetry=1, barrier=(21.0212,'kJ/mol'), semiclassical=False),
    ],
    spinMultiplicity = 6,
    opticalIsomers = 1,
    molecularWeight = (93.1033,'amu'),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[0.664631,0.0661338,-6.80905e-05,3.49188e-08,-6.9627e-12,110809,29.1364], Tmin=(100,'K'), Tmax=(1232.01,'K')), NASAPolynomial(coeffs=[16.5658,0.0145072,-5.23424e-06,9.06098e-10,-6.08448e-14,106891,-50.8953], Tmin=(1232.01,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(920.269,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(270.22,'J/(mol*K)'), comment="""Thermo group additivity estimation: group(O2s-(Cds-Cd)H) + group(Cs-(Cds-Cds)CsHH) + group(Cs-(Cds-Cds)CsHH) + group(Cds-CdsCsH) + group(Cds-CdsCsH) + group(Cds-CdsOsH) + group(Cds-CdsHH) + radical(Cds_S) + radical(C=COJ) + radical(Cds_S) + radical(CCJ2_triplet)"""),
)

species(
    label = '[CH2][CH]C=[C][C]=C[O](25022)',
    structure = SMILES('[CH2][CH]C=[C][C]=C[O]'),
    E0 = (712.27,'kJ/mol'),
    modes = [
        HarmonicOscillator(frequencies=([1670,1700,300,440,3025,407.5,1350,352.5,2995,3025,975,1000,1300,1375,400,500,1630,1680,3000,3100,440,815,1455,1000,200.078,200.265,200.914],'cm^-1')),
        HinderedRotor(inertia=(1.63594,'amu*angstrom^2'), symmetry=1, barrier=(46.5364,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(1.6166,'amu*angstrom^2'), symmetry=1, barrier=(46.5358,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(1.62583,'amu*angstrom^2'), symmetry=1, barrier=(46.5339,'kJ/mol'), semiclassical=False),
    ],
    spinMultiplicity = 6,
    opticalIsomers = 1,
    molecularWeight = (93.1033,'amu'),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[0.259313,0.0686097,-7.03502e-05,3.59266e-08,-6.93864e-12,85812.7,27.9761], Tmin=(100,'K'), Tmax=(1437.66,'K')), NASAPolynomial(coeffs=[17.6699,0.0116723,-2.07961e-06,1.57927e-10,-3.90246e-15,81684.6,-59.2866], Tmin=(1437.66,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(712.27,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(270.22,'J/(mol*K)'), comment="""Thermo group additivity estimation: group(O2s-(Cds-Cd)H) + group(Cs-(Cds-Cds)CsHH) + group(Cs-CsHHH) + group(Cds-CdsCsH) + group(Cds-Cds(Cds-Cds)H) + group(Cds-Cds(Cds-Cds)H) + group(Cds-CdsOsH) + radical(Allyl_S) + radical(C=CJC=C) + radical(C=CJC=C) + radical(RCCJ) + radical(C=COJ)"""),
)

species(
    label = 'C[C]C=[C][C]=C[O](25023)',
    structure = SMILES('C[C]C=[C][C]=C[O]'),
    E0 = (814.126,'kJ/mol'),
    modes = [
        HarmonicOscillator(frequencies=([2750,2800,2850,1350,1500,750,1050,1375,1000,2995,3025,975,1000,1300,1375,400,500,1630,1680,1670,1700,300,440,180,180,180,180],'cm^-1')),
        HinderedRotor(inertia=(1.64269,'amu*angstrom^2'), symmetry=1, barrier=(37.7687,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(1.64208,'amu*angstrom^2'), symmetry=1, barrier=(37.7547,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(1.6427,'amu*angstrom^2'), symmetry=1, barrier=(37.7689,'kJ/mol'), semiclassical=False),
    ],
    spinMultiplicity = 6,
    opticalIsomers = 1,
    molecularWeight = (93.1033,'amu'),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[0.104958,0.0746646,-8.35172e-05,4.60401e-08,-9.58186e-12,98066.5,27.057], Tmin=(100,'K'), Tmax=(1319.25,'K')), NASAPolynomial(coeffs=[18.5433,0.0107301,-1.69395e-06,7.84856e-11,2.17102e-15,93900.2,-64.3576], Tmin=(1319.25,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(814.126,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(270.22,'J/(mol*K)'), comment="""Thermo group additivity estimation: group(O2s-(Cds-Cd)H) + group(Cs-(Cds-Cds)CsHH) + group(Cs-CsHHH) + group(Cds-CdsCsH) + group(Cds-Cds(Cds-Cds)H) + group(Cds-Cds(Cds-Cds)H) + group(Cds-CdsOsH) + radical(CCJ2_triplet) + radical(C=CJC=C) + radical(C=CJC=C) + radical(C=COJ)"""),
)

species(
    label = 'C=C=CC=C=C[O](22623)',
    structure = SMILES('C=C=CC=C=C[O]'),
    E0 = (307.749,'kJ/mol'),
    modes = [
        HarmonicOscillator(frequencies=([2950,3100,1380,975,1025,1650,2995,3010,3025,975,987.5,1000,1300,1337.5,1375,400,450,500,1630,1655,1680,540,563.333,586.667,610,1970,2140,180,180],'cm^-1')),
        HinderedRotor(inertia=(1.44486,'amu*angstrom^2'), symmetry=1, barrier=(33.2202,'kJ/mol'), semiclassical=False),
    ],
    spinMultiplicity = 2,
    opticalIsomers = 1,
    molecularWeight = (93.1033,'amu'),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[0.822392,0.0563928,-2.37759e-05,-2.27834e-08,1.65251e-11,37140.4,22.3892], Tmin=(100,'K'), Tmax=(935.614,'K')), NASAPolynomial(coeffs=[20.0281,0.00835898,-1.398e-06,1.98487e-10,-1.72302e-14,32055.2,-76.9595], Tmin=(935.614,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(307.749,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(278.535,'J/(mol*K)'), comment="""Thermo group additivity estimation: group(O2s-(Cds-Cd)H) + group(Cds-Cds(Cds-Cds)H) + group(Cds-Cds(Cds-Cds)H) + group(Cds-CdsOsH) + group(Cds-CdsHH) + group(Cdd-CdsCds) + group(Cdd-CdsCds) + radical(C=COJ)"""),
)

species(
    label = 'C=C=CC1[C]C1[O](25024)',
    structure = SMILES('C=C=CC1[C]C1[O]'),
    E0 = (761.007,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
    molecularWeight = (93.1033,'amu'),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[1.37308,0.0434638,3.22731e-06,-4.23716e-08,2.11222e-11,91635.5,23.2975], Tmin=(100,'K'), Tmax=(958.357,'K')), NASAPolynomial(coeffs=[17.0156,0.01314,-4.03698e-06,7.51345e-10,-5.79993e-14,87031.6,-59.8803], Tmin=(958.357,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(761.007,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(278.535,'J/(mol*K)'), comment="""Thermo group additivity estimation: group(O2s-CsH) + group(Cs-(Cds-Cdd-Cd)CsCsH) + group(Cs-CsCsOsH) + group(Cs-CsCsHH) + group(Cds-CdsCsH) + group(Cds-CdsHH) + group(Cdd-CdsCds) + ring(Cyclopropane) + radical(CC(C)OJ) + radical(CCJ2_triplet)"""),
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
    E0 = (476.354,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS2',
    E0 = (598.397,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS3',
    E0 = (686.668,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS4',
    E0 = (692.465,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS5',
    E0 = (675.87,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS6',
    E0 = (665.835,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS7',
    E0 = (774.163,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS8',
    E0 = (690.818,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS9',
    E0 = (690.818,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS10',
    E0 = (669.787,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS11',
    E0 = (732.863,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS12',
    E0 = (688.617,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS13',
    E0 = (698.679,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS14',
    E0 = (698.679,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS15',
    E0 = (683.959,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS16',
    E0 = (721.787,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS17',
    E0 = (628.013,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS18',
    E0 = (566.927,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS19',
    E0 = (666.1,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS20',
    E0 = (887.154,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS21',
    E0 = (887.154,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS22',
    E0 = (867.847,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS23',
    E0 = (870.798,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS24',
    E0 = (614.003,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS25',
    E0 = (614.003,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS26',
    E0 = (552.638,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS27',
    E0 = (499.215,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS28',
    E0 = (569.448,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS29',
    E0 = (569.448,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS30',
    E0 = (501.327,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS31',
    E0 = (479.283,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS32',
    E0 = (484.638,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS33',
    E0 = (1084.62,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS34',
    E0 = (1043.27,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS35',
    E0 = (1115.73,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS36',
    E0 = (639.167,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS37',
    E0 = (652.894,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS38',
    E0 = (761.06,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS39',
    E0 = (693.871,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS40',
    E0 = (582.586,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS41',
    E0 = (884.525,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS42',
    E0 = (663.553,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS43',
    E0 = (614.409,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS44',
    E0 = (592.691,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS45',
    E0 = (545.997,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS46',
    E0 = (621.61,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS47',
    E0 = (501.327,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS48',
    E0 = (878.272,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS49',
    E0 = (789.489,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS50',
    E0 = (744.057,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS51',
    E0 = (777.702,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS52',
    E0 = (912.147,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS53',
    E0 = (815.368,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS54',
    E0 = (691.245,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS55',
    E0 = (881.125,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS56',
    E0 = (764.051,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS57',
    E0 = (501.327,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS58',
    E0 = (720.554,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS59',
    E0 = (768.928,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS60',
    E0 = (861.687,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS61',
    E0 = (539.571,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS62',
    E0 = (619.238,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS63',
    E0 = (737.163,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS64',
    E0 = (718.507,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS65',
    E0 = (910.116,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS66',
    E0 = (715.419,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS67',
    E0 = (994.046,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS68',
    E0 = (499.215,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS69',
    E0 = (479.283,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS70',
    E0 = (817.406,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS71',
    E0 = (809.391,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS72',
    E0 = (814.883,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS73',
    E0 = (948.122,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS74',
    E0 = (775.67,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS75',
    E0 = (839.099,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS76',
    E0 = (476.354,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS77',
    E0 = (761.007,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

reaction(
    label = 'reaction1',
    reactants = ['C=[C]C=C[C]=C[O](22614)'],
    products = ['C#CC=O(21959)', 'C#C[CH2](17441)'],
    transitionState = 'TS1',
    kinetics = Arrhenius(A=(5e+12,'s^-1'), n=0, Ea=(0,'kJ/mol'), T0=(1,'K'), Tmin=(300,'K'), Tmax=(1500,'K'), comment="""Exact match found for rate rule [RJJ]
Euclidian distance = 0
family: 1,4_Linear_birad_scission"""),
)

reaction(
    label = 'reaction2',
    reactants = ['C=[C]C=C[C]=C[O](22614)'],
    products = ['C=[C][CH]C1[C]=CO1(24975)'],
    transitionState = 'TS2',
    kinetics = Arrhenius(A=(2.724e+10,'s^-1'), n=0.478, Ea=(122.043,'kJ/mol'), T0=(1,'K'), comment="""Estimated using an average for rate rule [R5;doublebond_intra;radadd_intra_O]
Euclidian distance = 0
family: Intra_R_Add_Exocyclic"""),
)

reaction(
    label = 'reaction3',
    reactants = ['H(8)', '[CH2]C#CC=C=C[O](24976)'],
    products = ['C=[C]C=C[C]=C[O](22614)'],
    transitionState = 'TS3',
    kinetics = Arrhenius(A=(6.23e+08,'cm^3/(mol*s)'), n=1.429, Ea=(16.6816,'kJ/mol'), T0=(1,'K'), comment="""From training reaction 193 used for Ca_Ca;HJ
Exact match found for rate rule [Ca_Ca;HJ]
Euclidian distance = 0
family: R_Addition_MultipleBond"""),
)

reaction(
    label = 'reaction4',
    reactants = ['H(8)', '[CH2][C]=CC#CC=O(24977)'],
    products = ['C=[C]C=C[C]=C[O](22614)'],
    transitionState = 'TS4',
    kinetics = Arrhenius(A=(6.23e+08,'cm^3/(mol*s)'), n=1.429, Ea=(16.6816,'kJ/mol'), T0=(1,'K'), comment="""From training reaction 193 used for Ca_Ca;HJ
Exact match found for rate rule [Ca_Ca;HJ]
Euclidian distance = 0
family: R_Addition_MultipleBond"""),
)

reaction(
    label = 'reaction5',
    reactants = ['H(8)', 'C#CC=C[C]=C[O](22635)'],
    products = ['C=[C]C=C[C]=C[O](22614)'],
    transitionState = 'TS5',
    kinetics = Arrhenius(A=(2.156e+09,'cm^3/(mol*s)'), n=1.502, Ea=(9.92026,'kJ/mol'), T0=(1,'K'), comment="""From training reaction 194 used for Ct-H_Ct-Cd;HJ
Exact match found for rate rule [Ct-H_Ct-Cd;HJ]
Euclidian distance = 0
family: R_Addition_MultipleBond"""),
)

reaction(
    label = 'reaction6',
    reactants = ['H(8)', 'C=C=C[CH][C]=C=O(24978)'],
    products = ['C=[C]C=C[C]=C[O](22614)'],
    transitionState = 'TS6',
    kinetics = Arrhenius(A=(8.22e+08,'cm^3/(mol*s)'), n=1.533, Ea=(7.77387,'kJ/mol'), T0=(1,'K'), comment="""From training reaction 192 used for Cd_R;HJ
Exact match found for rate rule [Cd_R;HJ]
Euclidian distance = 0
family: R_Addition_MultipleBond"""),
)

reaction(
    label = 'reaction7',
    reactants = ['O(T)(63)', 'C#CC=C[C]=C(21245)'],
    products = ['C=[C]C=C[C]=C[O](22614)'],
    transitionState = 'TS7',
    kinetics = Arrhenius(A=(90.1275,'m^3/(mol*s)'), n=1.72743, Ea=(9.65468,'kJ/mol'), T0=(1,'K'), comment="""Estimated using template [Ct-H_Ct-Cd;YJ] for rate rule [Ct-H_Ct-Cd;O_atom_triplet]
Euclidian distance = 2.0
family: R_Addition_MultipleBond"""),
)

reaction(
    label = 'reaction8',
    reactants = ['C=[C]C=C[C]=C[O](22614)'],
    products = ['[CH2]C=[C]C=C=C[O](24979)'],
    transitionState = 'TS8',
    kinetics = Arrhenius(A=(1.806e+09,'s^-1'), n=1.172, Ea=(214.463,'kJ/mol'), T0=(1,'K'), comment="""Estimated using an average for rate rule [R2H_S;Cd_rad_out_Cd;Cd_H_out_doubleC]
Euclidian distance = 0
family: intra_H_migration"""),
)

reaction(
    label = 'reaction9',
    reactants = ['C=[C]C=[C]C=C[O](24980)'],
    products = ['C=[C]C=C[C]=C[O](22614)'],
    transitionState = 'TS9',
    kinetics = Arrhenius(A=(1.806e+09,'s^-1'), n=1.172, Ea=(214.463,'kJ/mol'), T0=(1,'K'), comment="""Estimated using an average for rate rule [R2H_S;Cd_rad_out_Cd;Cd_H_out_doubleC]
Euclidian distance = 0
family: intra_H_migration"""),
)

reaction(
    label = 'reaction10',
    reactants = ['[CH]C=CC=C=C[O](24787)'],
    products = ['C=[C]C=C[C]=C[O](22614)'],
    transitionState = 'TS10',
    kinetics = Arrhenius(A=(481900,'s^-1'), n=2.375, Ea=(167.958,'kJ/mol'), T0=(1,'K'), comment="""From training reaction 121 used for R2H_D;Cd_rad_out_singleH;Cd_H_out_singleDe
Exact match found for rate rule [R2H_D;Cd_rad_out_singleH;Cd_H_out_singleDe]
Euclidian distance = 0
family: intra_H_migration"""),
)

reaction(
    label = 'reaction11',
    reactants = ['C=[C]C=C[C]=[C]O(24981)'],
    products = ['C=[C]C=C[C]=C[O](22614)'],
    transitionState = 'TS11',
    kinetics = Arrhenius(A=(6.25466e+06,'s^-1'), n=1.80084, Ea=(158.227,'kJ/mol'), T0=(1,'K'), comment="""Estimated using average of templates [R2H_S;Y_rad_out;O_H_out] + [R2H_S;Cd_rad_out;XH_out] for rate rule [R2H_S;Cd_rad_out;O_H_out]
Euclidian distance = 1.0
family: intra_H_migration"""),
)

reaction(
    label = 'reaction12',
    reactants = ['C=[C]C=C[C]=C[O](22614)'],
    products = ['[CH2][C]=CC=C[C]=O(24982)'],
    transitionState = 'TS12',
    kinetics = Arrhenius(A=(1.0945e+07,'s^-1'), n=1.951, Ea=(212.263,'kJ/mol'), T0=(1,'K'), comment="""Estimated using an average for rate rule [R2H_D;Cd_rad_out_singleDe_Cd;XH_out]
Euclidian distance = 0
family: intra_H_migration"""),
)

reaction(
    label = 'reaction13',
    reactants = ['C=CC=[C][C]=C[O](24983)'],
    products = ['C=[C]C=C[C]=C[O](22614)'],
    transitionState = 'TS13',
    kinetics = Arrhenius(A=(191.5,'s^-1'), n=3.05, Ea=(222.325,'kJ/mol'), T0=(1,'K'), comment="""Estimated using an average for rate rule [R3H_DS;Cd_rad_out;Cd_H_out_doubleC]
Euclidian distance = 0
family: intra_H_migration"""),
)

reaction(
    label = 'reaction14',
    reactants = ['C=[C][C]=CC=C[O](24984)'],
    products = ['C=[C]C=C[C]=C[O](22614)'],
    transitionState = 'TS14',
    kinetics = Arrhenius(A=(191.5,'s^-1'), n=3.05, Ea=(222.325,'kJ/mol'), T0=(1,'K'), comment="""Estimated using an average for rate rule [R3H_DS;Cd_rad_out;Cd_H_out_doubleC]
Euclidian distance = 0
family: intra_H_migration"""),
)

reaction(
    label = 'reaction15',
    reactants = ['C=[C]C=[C][C]=CO(24985)'],
    products = ['C=[C]C=C[C]=C[O](22614)'],
    transitionState = 'TS15',
    kinetics = Arrhenius(A=(1.13341e+09,'s^-1'), n=1.04717, Ea=(150.072,'kJ/mol'), T0=(1,'K'), comment="""Estimated using template [R4Hall;Cd_rad_out_Cd;XH_out] for rate rule [R4HJ_1;Cd_rad_out_Cd;O_H_out]
Euclidian distance = 1.41421356237
family: intra_H_migration"""),
)

reaction(
    label = 'reaction16',
    reactants = ['C=[C]C=C[C]=C[O](22614)'],
    products = ['[CH]=C=C[CH]C=C[O](21998)'],
    transitionState = 'TS16',
    kinetics = Arrhenius(A=(2.04e+11,'s^-1'), n=0.627, Ea=(245.433,'kJ/mol'), T0=(1,'K'), comment="""Estimated using template [RnH;Cd_rad_out_Cd;Cd_H_out_singleH] for rate rule [R5HJ_3;Cd_rad_out_Cd;Cd_H_out_singleH]
Euclidian distance = 2.0
Multiplied by reaction path degeneracy 2.0
family: intra_H_migration"""),
)

reaction(
    label = 'reaction17',
    reactants = ['C=[C]C=C[C]=C[O](22614)'],
    products = ['[CH2]C=C[CH][C]=C=O(24986)'],
    transitionState = 'TS17',
    kinetics = Arrhenius(A=(4.81182e+08,'s^-1'), n=1.25566, Ea=(151.659,'kJ/mol'), T0=(1,'K'), comment="""Estimated using template [RnH;Cd_rad_out_Cd;XH_out] for rate rule [R5HJ_3;Cd_rad_out_Cd;XH_out]
Euclidian distance = 2.0
family: intra_H_migration"""),
)

reaction(
    label = 'reaction18',
    reactants = ['C=[C][C]=C[C]=CO(24987)'],
    products = ['C=[C]C=C[C]=C[O](22614)'],
    transitionState = 'TS18',
    kinetics = Arrhenius(A=(136000,'s^-1'), n=1.9199, Ea=(33.0402,'kJ/mol'), T0=(1,'K'), comment="""Estimated using template [R5Hall;Cd_rad_out;XH_out] for rate rule [R5HJ_2;Cd_rad_out;O_H_out]
Euclidian distance = 1.41421356237
family: intra_H_migration"""),
)

reaction(
    label = 'reaction19',
    reactants = ['[CH]=C=C[CH][C]=CO(24792)'],
    products = ['C=[C]C=C[C]=C[O](22614)'],
    transitionState = 'TS19',
    kinetics = Arrhenius(A=(1.86943e+06,'s^-1'), n=1.85754, Ea=(151.521,'kJ/mol'), T0=(1,'K'), comment="""Estimated using template [RnH;Cd_rad_out_singleH;XH_out] for rate rule [R7Hall;Cd_rad_out_singleH;O_H_out]
Euclidian distance = 1.41421356237
family: intra_H_migration"""),
)

reaction(
    label = 'reaction20',
    reactants = ['H(8)', 'C=[C][C]=C[C]=C[O](24988)'],
    products = ['C=[C]C=C[C]=C[O](22614)'],
    transitionState = 'TS20',
    kinetics = Arrhenius(A=(4.34078e+06,'m^3/(mol*s)'), n=0.278577, Ea=(0,'kJ/mol'), T0=(1,'K'), comment="""Estimated using an average for rate rule [Y_rad;H_rad]
Euclidian distance = 0
family: R_Recombination
Ea raised from -1.4 to 0 kJ/mol."""),
)

reaction(
    label = 'reaction21',
    reactants = ['H(8)', 'C=[C]C=[C][C]=C[O](24989)'],
    products = ['C=[C]C=C[C]=C[O](22614)'],
    transitionState = 'TS21',
    kinetics = Arrhenius(A=(4.34078e+06,'m^3/(mol*s)'), n=0.278577, Ea=(0,'kJ/mol'), T0=(1,'K'), comment="""Estimated using an average for rate rule [Y_rad;H_rad]
Euclidian distance = 0
family: R_Recombination
Ea raised from -1.4 to 0 kJ/mol."""),
)

reaction(
    label = 'reaction22',
    reactants = ['H(8)', '[CH]=C=C[CH][C]=C[O](22638)'],
    products = ['C=[C]C=C[C]=C[O](22614)'],
    transitionState = 'TS22',
    kinetics = Arrhenius(A=(4.34078e+06,'m^3/(mol*s)'), n=0.278577, Ea=(0,'kJ/mol'), T0=(1,'K'), comment="""Estimated using an average for rate rule [Y_rad;H_rad]
Euclidian distance = 0
family: R_Recombination
Ea raised from -1.4 to 0 kJ/mol."""),
)

reaction(
    label = 'reaction23',
    reactants = ['H(8)', '[CH2][C]=C[CH][C]=C=O(24990)'],
    products = ['C=[C]C=C[C]=C[O](22614)'],
    transitionState = 'TS23',
    kinetics = Arrhenius(A=(4.34078e+06,'m^3/(mol*s)'), n=0.278577, Ea=(0,'kJ/mol'), T0=(1,'K'), comment="""Estimated using an average for rate rule [H_rad;Y_rad]
Euclidian distance = 0
family: R_Recombination
Ea raised from -1.4 to 0 kJ/mol."""),
)

reaction(
    label = 'reaction24',
    reactants = ['C=[C]C=C[C]=C[O](22614)'],
    products = ['C=C1[CH]C1[C]=C[O](24883)'],
    transitionState = 'TS24',
    kinetics = Arrhenius(A=(5.605e+12,'s^-1'), n=0.275, Ea=(137.649,'kJ/mol'), T0=(1,'K'), comment="""Estimated using an average for rate rule [R3_D;doublebond_intra_pri;radadd_intra_cddouble]
Euclidian distance = 0
family: Intra_R_Add_Endocyclic"""),
)

reaction(
    label = 'reaction25',
    reactants = ['C=[C]C=C[C]=C[O](22614)'],
    products = ['C=[C]C1[CH]C1=C[O](24969)'],
    transitionState = 'TS25',
    kinetics = Arrhenius(A=(5.605e+12,'s^-1'), n=0.275, Ea=(137.649,'kJ/mol'), T0=(1,'K'), comment="""Estimated using an average for rate rule [R3_D;doublebond_intra_pri;radadd_intra_cddouble]
Euclidian distance = 0
family: Intra_R_Add_Endocyclic"""),
)

reaction(
    label = 'reaction26',
    reactants = ['C=[C]C=C[C]=C[O](22614)'],
    products = ['C=[C]C1[CH][C]=CO1(24991)'],
    transitionState = 'TS26',
    kinetics = Arrhenius(A=(1.72761e+09,'s^-1'), n=0.568593, Ea=(76.2838,'kJ/mol'), T0=(1,'K'), comment="""Estimated using template [R5_linear;doublebond_intra_pri;radadd_intra] for rate rule [R5_linear;doublebond_intra_pri;radadd_intra_O]
Euclidian distance = 1.0
family: Intra_R_Add_Endocyclic"""),
)

reaction(
    label = 'reaction27',
    reactants = ['C=[C]C=C[C]=C[O](22614)'],
    products = ['C=C=C[CH]C=C=O(24992)'],
    transitionState = 'TS27',
    kinetics = Arrhenius(A=(1.949e+11,'s^-1'), n=0.486, Ea=(22.8614,'kJ/mol'), T0=(1,'K'), comment="""Estimated using an average for rate rule [R2radExo;Y_rad_De;XH_Rrad]
Euclidian distance = 0
family: Intra_Disproportionation"""),
)

reaction(
    label = 'reaction28',
    reactants = ['C=[C]C=C[C]=C[O](22614)'],
    products = ['[CH2]C#CC=CC=O(24993)'],
    transitionState = 'TS28',
    kinetics = Arrhenius(A=(5.4e+09,'s^-1'), n=-0.305, Ea=(93.094,'kJ/mol'), T0=(1,'K'), Tmin=(300,'K'), Tmax=(1500,'K'), comment="""Estimated using template [R3;Y_rad_De;XH_Rrad_De] for rate rule [R3radExo;Y_rad_De;XH_Rrad_De]
Euclidian distance = 1.0
family: Intra_Disproportionation"""),
)

reaction(
    label = 'reaction29',
    reactants = ['C=[C]C=C[C]=C[O](22614)'],
    products = ['[CH2]C=CC#CC=O(24994)'],
    transitionState = 'TS29',
    kinetics = Arrhenius(A=(5.4e+09,'s^-1'), n=-0.305, Ea=(93.094,'kJ/mol'), T0=(1,'K'), Tmin=(300,'K'), Tmax=(1500,'K'), comment="""Estimated using template [R3;Y_rad_De;XH_Rrad_De] for rate rule [R3radExo;Y_rad_De;XH_Rrad_De]
Euclidian distance = 1.0
family: Intra_Disproportionation"""),
)

reaction(
    label = 'reaction30',
    reactants = ['C=[C]C=C[C]=C[O](22614)'],
    products = ['C=[C]C=C=C=CO(24995)'],
    transitionState = 'TS30',
    kinetics = Arrhenius(A=(5.14222e+08,'s^-1'), n=0.311, Ea=(24.9733,'kJ/mol'), T0=(1,'K'), comment="""Estimated using template [R4;Y_rad;XH_Rrad] for rate rule [R4radEndo;Y_rad;XH_Rrad]
Euclidian distance = 1.0
family: Intra_Disproportionation"""),
)

reaction(
    label = 'reaction31',
    reactants = ['C=[C]C=C[C]=C[O](22614)'],
    products = ['C=C=CC=C1[CH]O1(24996)'],
    transitionState = 'TS31',
    kinetics = Arrhenius(A=(2.18e+16,'s^-1'), n=0, Ea=(2.9288,'kJ/mol'), T0=(1,'K'), comment="""Estimated using template [R3;Y_rad_out;Ypri_rad_out] for rate rule [R3_SD;O_rad;CdsinglepriDe_rad_out]
Euclidian distance = 2.44948974278
family: Birad_recombination"""),
)

reaction(
    label = 'reaction32',
    reactants = ['C=[C]C=C[C]=C[O](22614)'],
    products = ['C=C1C=C[C]1C=O(24864)'],
    transitionState = 'TS32',
    kinetics = Arrhenius(A=(1.62e+12,'s^-1'), n=-0.305, Ea=(8.28432,'kJ/mol'), T0=(1,'K'), comment="""Estimated using template [R4;Y_rad_out;Ypri_rad_out] for rate rule [R4_SDS;Y_rad_out;Ypri_rad_out]
Euclidian distance = 1.0
family: Birad_recombination"""),
)

reaction(
    label = 'reaction33',
    reactants = ['[C]=C(584)', '[CH]C=C=C[O](24997)'],
    products = ['C=[C]C=C[C]=C[O](22614)'],
    transitionState = 'TS33',
    kinetics = Arrhenius(A=(1.14854e+06,'m^3/(mol*s)'), n=0.575199, Ea=(34.3157,'kJ/mol'), T0=(1,'K'), comment="""Estimated using template [Y_rad;Birad] for rate rule [Cd_pri_rad;Birad]
Euclidian distance = 2.0
family: Birad_R_Recombination"""),
)

reaction(
    label = 'reaction34',
    reactants = ['[CH]C=C=C(5316)', '[C]C=O(24998)'],
    products = ['C=[C]C=C[C]=C[O](22614)'],
    transitionState = 'TS34',
    kinetics = Arrhenius(A=(1.14854e+06,'m^3/(mol*s)'), n=0.575199, Ea=(34.3157,'kJ/mol'), T0=(1,'K'), comment="""Estimated using template [Y_rad;Birad] for rate rule [Cd_pri_rad;Birad]
Euclidian distance = 2.0
family: Birad_R_Recombination"""),
)

reaction(
    label = 'reaction35',
    reactants = ['O(T)(63)', '[CH]=C=CC=[C][CH2](19899)'],
    products = ['C=[C]C=C[C]=C[O](22614)'],
    transitionState = 'TS35',
    kinetics = Arrhenius(A=(93609.6,'m^3/(mol*s)'), n=1.13083, Ea=(163.847,'kJ/mol'), T0=(1,'K'), Tmin=(303.03,'K'), Tmax=(2000,'K'), comment="""From training reaction 2 used for Y_rad;O_birad
Exact match found for rate rule [Y_rad;O_birad]
Euclidian distance = 0
family: Birad_R_Recombination"""),
)

reaction(
    label = 'reaction36',
    reactants = ['C=[C]C=C[C]=C[O](22614)'],
    products = ['[O]C=[C]C1C=[C]C1(24785)'],
    transitionState = 'TS36',
    kinetics = Arrhenius(A=(4.713e+10,'s^-1'), n=0.481, Ea=(162.813,'kJ/mol'), T0=(1,'K'), comment="""Estimated using an average for rate rule [R5;doublebond_intra;radadd_intra_cs2H]
Euclidian distance = 0
family: Intra_R_Add_Exocyclic
Ea raised from 160.1 to 162.8 kJ/mol to match endothermicity of reaction."""),
)

reaction(
    label = 'reaction37',
    reactants = ['[CH]=C=C[O](8556)', 'C#C[CH2](17441)'],
    products = ['C=[C]C=C[C]=C[O](22614)'],
    transitionState = 'TS37',
    kinetics = Arrhenius(A=(0.74854,'m^3/(mol*s)'), n=1.927, Ea=(54.5384,'kJ/mol'), T0=(1,'K'), comment="""Estimated using an average for rate rule [Ct_Ct;CdsJ=Cdd]
Euclidian distance = 0
family: R_Addition_MultipleBond"""),
)

reaction(
    label = 'reaction38',
    reactants = ['C[C]=[C]C=C=C[O](24999)'],
    products = ['C=[C]C=C[C]=C[O](22614)'],
    transitionState = 'TS38',
    kinetics = Arrhenius(A=(3.85113e+09,'s^-1'), n=1.0541, Ea=(193.078,'kJ/mol'), T0=(1,'K'), Tmin=(300,'K'), Tmax=(1600,'K'), comment="""Estimated using template [R3Hall;Cd_rad_out_singleDe_Cd;Cs_H_out_2H] for rate rule [R3HJ;Cd_rad_out_singleDe_Cd;Cs_H_out_2H]
Euclidian distance = 1.0
Multiplied by reaction path degeneracy 3.0
family: intra_H_migration"""),
)

reaction(
    label = 'reaction39',
    reactants = ['C[C][CH]C#CC=O(25000)'],
    products = ['C=[C]C=C[C]=C[O](22614)'],
    transitionState = 'TS39',
    kinetics = Arrhenius(A=(3.33e+08,'s^-1'), n=1.1915, Ea=(103.605,'kJ/mol'), T0=(1,'K'), comment="""Estimated using template [R4Hall;Cd_rad_out_double;Cs_H_out_2H] for rate rule [R4HJ_2;Cd_rad_out_double;Cs_H_out_2H]
Euclidian distance = 1.0
Multiplied by reaction path degeneracy 3.0
family: intra_H_migration"""),
)

reaction(
    label = 'reaction40',
    reactants = ['C[C]=C[CH][C]=C=O(25001)'],
    products = ['C=[C]C=C[C]=C[O](22614)'],
    transitionState = 'TS40',
    kinetics = Arrhenius(A=(3.31883e+06,'s^-1'), n=1.02765, Ea=(75.0925,'kJ/mol'), T0=(1,'K'), comment="""Estimated using template [R6Hall;Y_rad_out;Cs_H_out_2H] for rate rule [R6HJ_4;Cd_rad_out;Cs_H_out_2H]
Euclidian distance = 1.41421356237
Multiplied by reaction path degeneracy 3.0
family: intra_H_migration"""),
)

reaction(
    label = 'reaction41',
    reactants = ['[CH]=C=C[O](8556)', '[CH][C]=C(18825)'],
    products = ['C=[C]C=C[C]=C[O](22614)'],
    transitionState = 'TS41',
    kinetics = Arrhenius(A=(40577.9,'m^3/(mol*s)'), n=0.702818, Ea=(0,'kJ/mol'), T0=(1,'K'), comment="""Estimated using an average for rate rule [Cd_allenic;Y_rad]
Euclidian distance = 0
family: R_Recombination
Ea raised from -11.4 to 0 kJ/mol."""),
)

reaction(
    label = 'reaction42',
    reactants = ['C=[C]C=C[C]=C[O](22614)'],
    products = ['[CH2][C]=CC=C1[CH]O1(25002)'],
    transitionState = 'TS42',
    kinetics = Arrhenius(A=(6.59543e+11,'s^-1'), n=0.400725, Ea=(187.199,'kJ/mol'), T0=(1,'K'), comment="""Estimated using template [R3_D;doublebond_intra;radadd_intra] for rate rule [R3_D;doublebond_intra;radadd_intra_O]
Euclidian distance = 1.0
family: Intra_R_Add_Endocyclic"""),
)

reaction(
    label = 'reaction43',
    reactants = ['C=[C]C=C[C]=C[O](22614)'],
    products = ['[CH2]C1=C[CH]C1=C[O](25003)'],
    transitionState = 'TS43',
    kinetics = Arrhenius(A=(1.953e+11,'s^-1'), n=0.387, Ea=(138.055,'kJ/mol'), T0=(1,'K'), comment="""Estimated using an average for rate rule [R4_D_D;doublebond_intra;radadd_intra]
Euclidian distance = 0
family: Intra_R_Add_Endocyclic"""),
)

reaction(
    label = 'reaction44',
    reactants = ['C=[C]C=C[C]=C[O](22614)'],
    products = ['[O]C=C1[CH]C=[C]C1(25004)'],
    transitionState = 'TS44',
    kinetics = Arrhenius(A=(4.4946e+10,'s^-1'), n=0.314866, Ea=(116.337,'kJ/mol'), T0=(1,'K'), comment="""Estimated using an average for rate rule [R5_linear;doublebond_intra;radadd_intra_cs2H]
Euclidian distance = 0
family: Intra_R_Add_Endocyclic"""),
)

reaction(
    label = 'reaction45',
    reactants = ['C=[C]C=C[C]=C[O](22614)'],
    products = ['[CH2]C1=CC=[C]C1[O](24930)'],
    transitionState = 'TS45',
    kinetics = Arrhenius(A=(6.01304e+12,'s^-1'), n=-0.3725, Ea=(69.6427,'kJ/mol'), T0=(1,'K'), comment="""Estimated using an average for rate rule [R5_MS;doublebond_intra_CdCdd;radadd_intra]
Euclidian distance = 0
family: Intra_R_Add_Endocyclic"""),
)

reaction(
    label = 'reaction46',
    reactants = ['C=[C]C=C[C]=C[O](22614)'],
    products = ['[O]C1[C]=CC=[C]C1(24844)'],
    transitionState = 'TS46',
    kinetics = Arrhenius(A=(4.26053e+10,'s^-1'), n=0.2505, Ea=(145.256,'kJ/mol'), T0=(1,'K'), comment="""Estimated using an average for rate rule [R6_linear;doublebond_intra_CdCdd;radadd_intra_cs2H]
Euclidian distance = 0
family: Intra_R_Add_Endocyclic
Ea raised from 140.6 to 145.3 kJ/mol to match endothermicity of reaction."""),
)

reaction(
    label = 'reaction47',
    reactants = ['C=[C]C=C[C]=C[O](22614)'],
    products = ['[CH2]C=CC=C=C=O(25005)'],
    transitionState = 'TS47',
    kinetics = Arrhenius(A=(2.1261e+09,'s^-1'), n=0.137, Ea=(24.9733,'kJ/mol'), T0=(1,'K'), comment="""Estimated using template [R5;Y_rad;XH_Rrad] for rate rule [R5radExo;Y_rad;XH_Rrad]
Euclidian distance = 1.0
family: Intra_Disproportionation"""),
)

reaction(
    label = 'reaction48',
    reactants = ['[CH2][C]=[C]C[C]=C[O](25006)'],
    products = ['C=[C]C=C[C]=C[O](22614)'],
    transitionState = 'TS48',
    kinetics = Arrhenius(A=(2.94659e+10,'s^-1'), n=0.2847, Ea=(27.8529,'kJ/mol'), T0=(1,'K'), comment="""Estimated using average of templates [Rn;Y_rad;XH_Rrad_De] + [R2radExo;Y_rad;XH_Rrad] for rate rule [R2radExo;Y_rad;XH_Rrad_De]
Euclidian distance = 1.0
Multiplied by reaction path degeneracy 2.0
family: Intra_Disproportionation"""),
)

reaction(
    label = 'reaction49',
    reactants = ['[CH2][C]=C[C]=C[CH][O](25007)'],
    products = ['C=[C]C=C[C]=C[O](22614)'],
    transitionState = 'TS49',
    kinetics = Arrhenius(A=(1.949e+11,'s^-1'), n=0.486, Ea=(22.8614,'kJ/mol'), T0=(1,'K'), comment="""Estimated using an average for rate rule [R2radExo;Y_rad_De;XH_Rrad]
Euclidian distance = 0
family: Intra_Disproportionation"""),
)

reaction(
    label = 'reaction50',
    reactants = ['[CH2][C]=C[CH]C=[C][O](25008)'],
    products = ['C=[C]C=C[C]=C[O](22614)'],
    transitionState = 'TS50',
    kinetics = Arrhenius(A=(1.4733e+10,'s^-1'), n=0.2847, Ea=(27.8529,'kJ/mol'), T0=(1,'K'), comment="""Estimated using average of templates [Rn;Y_rad;XH_Rrad_De] + [R2radExo;Y_rad;XH_Rrad] for rate rule [R2radExo;Y_rad;XH_Rrad_De]
Euclidian distance = 1.0
family: Intra_Disproportionation"""),
)

reaction(
    label = 'reaction51',
    reactants = ['[CH2][C]=[C][CH]C=C[O](25009)'],
    products = ['C=[C]C=C[C]=C[O](22614)'],
    transitionState = 'TS51',
    kinetics = Arrhenius(A=(7.437e+08,'s^-1'), n=1.045, Ea=(63.4002,'kJ/mol'), T0=(1,'K'), comment="""Estimated using an average for rate rule [R3radExo;Y_rad;XH_Rrad]
Euclidian distance = 0
family: Intra_Disproportionation"""),
)

reaction(
    label = 'reaction52',
    reactants = ['[CH2][C]=[C]C=[C]C[O](25010)'],
    products = ['C=[C]C=C[C]=C[O](22614)'],
    transitionState = 'TS52',
    kinetics = Arrhenius(A=(1.02844e+09,'s^-1'), n=0.311, Ea=(24.9733,'kJ/mol'), T0=(1,'K'), comment="""Estimated using template [R4;Y_rad;XH_Rrad] for rate rule [R4radEndo;Y_rad;XH_Rrad]
Euclidian distance = 1.0
Multiplied by reaction path degeneracy 2.0
family: Intra_Disproportionation"""),
)

reaction(
    label = 'reaction53',
    reactants = ['C=[C]C[C]=C=C[O](24306)'],
    products = ['C=[C]C=C[C]=C[O](22614)'],
    transitionState = 'TS53',
    kinetics = Arrhenius(A=(7.80481e+09,'s^-1'), n=1.09833, Ea=(177.686,'kJ/mol'), T0=(1,'K'), comment="""Estimated using an average for rate rule [R2H_S;Cd_rad_out_double;XH_out]
Euclidian distance = 0
Multiplied by reaction path degeneracy 2.0
family: intra_H_migration"""),
)

reaction(
    label = 'reaction10',
    reactants = ['[CH]=[C]CC=C=C[O](22626)'],
    products = ['C=[C]C=C[C]=C[O](22614)'],
    transitionState = 'TS54',
    kinetics = Arrhenius(A=(74200,'s^-1'), n=2.23, Ea=(44.3086,'kJ/mol'), T0=(1,'K'), comment="""Estimated using template [RnH;Cd_rad_out_singleH;Cs_H_out_H/Cd] for rate rule [R3HJ;Cd_rad_out_singleH;Cs_H_out_H/Cd]
Euclidian distance = 2.0
Multiplied by reaction path degeneracy 2.0
family: intra_H_migration"""),
)

reaction(
    label = 'reaction55',
    reactants = ['C=[C]C[CH][C]=C=O(24307)'],
    products = ['C=[C]C=C[C]=C[O](22614)'],
    transitionState = 'TS55',
    kinetics = Arrhenius(A=(6586.33,'s^-1'), n=2.5965, Ea=(315.162,'kJ/mol'), T0=(1,'K'), comment="""Estimated using an average for rate rule [R4H_MMS;Cd_rad_out;XH_out]
Euclidian distance = 0
Multiplied by reaction path degeneracy 2.0
family: intra_H_migration"""),
)

reaction(
    label = 'reaction56',
    reactants = ['C=[C]C=C[C]=C[O](22614)'],
    products = ['C=[C]C1C=[C]C1[O](25011)'],
    transitionState = 'TS56',
    kinetics = Arrhenius(A=(3.01156e+10,'s^-1'), n=0.428741, Ea=(287.697,'kJ/mol'), T0=(1,'K'), comment="""Estimated using average of templates [Rn;doublebond_intra_CdCdd;radadd_intra_cs] + [R4;doublebond_intra;radadd_intra_cs] for rate rule [R4;doublebond_intra_CdCdd;radadd_intra_cs]
Euclidian distance = 1.0
family: Intra_R_Add_Endocyclic
Ea raised from 285.0 to 287.7 kJ/mol to match endothermicity of reaction."""),
)

reaction(
    label = 'reaction57',
    reactants = ['C=[C]C=C[C]=C[O](22614)'],
    products = ['C=[C]CC=C=C=O(24300)'],
    transitionState = 'TS57',
    kinetics = Arrhenius(A=(5.14222e+08,'s^-1'), n=0.311, Ea=(24.9733,'kJ/mol'), T0=(1,'K'), comment="""Estimated using template [R4;Y_rad;XH_Rrad] for rate rule [R4radExo;Y_rad;XH_Rrad]
Euclidian distance = 1.0
family: Intra_Disproportionation"""),
)

reaction(
    label = 'reaction58',
    reactants = ['[CH2]C#CC[C]=C[O](25012)'],
    products = ['C=[C]C=C[C]=C[O](22614)'],
    transitionState = 'TS58',
    kinetics = Arrhenius(A=(7.80481e+09,'s^-1'), n=1.09833, Ea=(177.686,'kJ/mol'), T0=(1,'K'), comment="""Estimated using an average for rate rule [R2H_S;Cd_rad_out_double;XH_out]
Euclidian distance = 0
Multiplied by reaction path degeneracy 2.0
family: intra_H_migration"""),
)

reaction(
    label = 'reaction59',
    reactants = ['C=C=CC[C][C]=O(25013)'],
    products = ['C=[C]C=C[C]=C[O](22614)'],
    transitionState = 'TS59',
    kinetics = Arrhenius(A=(1.34586e+06,'s^-1'), n=1.99734, Ea=(143.275,'kJ/mol'), T0=(1,'K'), comment="""Estimated using template [R3Hall;Y_rad_out;Cs_H_out_H/Cd] for rate rule [R3HJ;Cd_rad_out;Cs_H_out_H/Cd]
Euclidian distance = 1.41421356237
Multiplied by reaction path degeneracy 2.0
family: intra_H_migration"""),
)

reaction(
    label = 'reaction60',
    reactants = ['C#C[CH]C[C]=C[O](25014)'],
    products = ['C=[C]C=C[C]=C[O](22614)'],
    transitionState = 'TS60',
    kinetics = Arrhenius(A=(6586.33,'s^-1'), n=2.5965, Ea=(315.162,'kJ/mol'), T0=(1,'K'), comment="""Estimated using an average for rate rule [R4H_MMS;Cd_rad_out_singleH;XH_out]
Euclidian distance = 0
Multiplied by reaction path degeneracy 2.0
family: intra_H_migration"""),
)

reaction(
    label = 'reaction61',
    reactants = ['C=[C]C=C[C]=C[O](22614)'],
    products = ['C=C1[CH]C=[C][CH]O1(24914)'],
    transitionState = 'TS61',
    kinetics = Arrhenius(A=(2.49515e+10,'s^-1'), n=0.243684, Ea=(63.2167,'kJ/mol'), T0=(1,'K'), comment="""Estimated using template [R6_linear;doublebond_intra;radadd_intra] for rate rule [R6_linear;doublebond_intra;radadd_intra_O]
Euclidian distance = 1.0
family: Intra_R_Add_Endocyclic"""),
)

reaction(
    label = 'reaction62',
    reactants = ['C=[C]C=C[C]=C[O](22614)'],
    products = ['[C]1=C[CH][C]=COC1(24829)'],
    transitionState = 'TS62',
    kinetics = Arrhenius(A=(3.22e+12,'s^-1'), n=-0.622, Ea=(142.884,'kJ/mol'), T0=(1,'K'), comment="""Estimated using template [R7;doublebond_intra_CdCdd;radadd_intra] for rate rule [R7;doublebond_intra_CdCdd;radadd_intra_O]
Euclidian distance = 1.0
family: Intra_R_Add_Endocyclic"""),
)

reaction(
    label = 'reaction63',
    reactants = ['[CH2]C=[C][CH][C]=C[O](25015)'],
    products = ['C=[C]C=C[C]=C[O](22614)'],
    transitionState = 'TS63',
    kinetics = Arrhenius(A=(1.949e+11,'s^-1'), n=0.486, Ea=(22.8614,'kJ/mol'), T0=(1,'K'), comment="""Estimated using an average for rate rule [R2radExo;Y_rad;XH_Rrad]
Euclidian distance = 0
family: Intra_Disproportionation"""),
)

reaction(
    label = 'reaction64',
    reactants = ['[CH]C=C[CH][C]=C[O](25016)'],
    products = ['C=[C]C=C[C]=C[O](22614)'],
    transitionState = 'TS64',
    kinetics = Arrhenius(A=(1.949e+11,'s^-1'), n=0.486, Ea=(22.8614,'kJ/mol'), T0=(1,'K'), comment="""Estimated using an average for rate rule [R2radExo;Y_rad;XH_Rrad]
Euclidian distance = 0
family: Intra_Disproportionation"""),
)

reaction(
    label = 'reaction65',
    reactants = ['[CH][C]=C[O](21209)', 'C#C[CH2](17441)'],
    products = ['C=[C]C=C[C]=C[O](22614)'],
    transitionState = 'TS65',
    kinetics = Arrhenius(A=(1.14854e+06,'m^3/(mol*s)'), n=0.575199, Ea=(34.3157,'kJ/mol'), T0=(1,'K'), comment="""Estimated using template [Y_rad;Birad] for rate rule [Cd_pri_rad;Birad]
Euclidian distance = 2.0
family: Birad_R_Recombination"""),
)

reaction(
    label = 'reaction66',
    reactants = ['C#CC=O(21959)', '[CH][C]=C(18825)'],
    products = ['C=[C]C=C[C]=C[O](22614)'],
    transitionState = 'TS66',
    kinetics = Arrhenius(A=(0.106247,'m^3/(mol*s)'), n=2.32278, Ea=(16.475,'kJ/mol'), T0=(1,'K'), comment="""Estimated using template [Ct-H_Ct-De;CJ] for rate rule [Ct-H_Ct-CO;CJ]
Euclidian distance = 1.0
family: R_Addition_MultipleBond"""),
)

reaction(
    label = 'reaction67',
    reactants = ['[CH]=O(373)', '[C]C=C[C]=C(19282)'],
    products = ['C=[C]C=C[C]=C[O](22614)'],
    transitionState = 'TS67',
    kinetics = Arrhenius(A=(1.14854e+06,'m^3/(mol*s)'), n=0.575199, Ea=(34.3157,'kJ/mol'), T0=(1,'K'), comment="""Estimated using template [Y_rad;Birad] for rate rule [CO_pri_rad;Birad]
Euclidian distance = 2.0
family: Birad_R_Recombination"""),
)

reaction(
    label = 'reaction68',
    reactants = ['C=[C]C=C[C]=C[O](22614)'],
    products = ['C=[C]C=C=CC=O(25017)'],
    transitionState = 'TS68',
    kinetics = Arrhenius(A=(1.949e+11,'s^-1'), n=0.486, Ea=(22.8614,'kJ/mol'), T0=(1,'K'), comment="""Estimated using an average for rate rule [R2radExo;Y_rad_De;XH_Rrad]
Euclidian distance = 0
family: Intra_Disproportionation"""),
)

reaction(
    label = 'reaction69',
    reactants = ['C=[C]C=C[C]=C[O](22614)'],
    products = ['C=[C]C1C=C1C=O(24974)'],
    transitionState = 'TS69',
    kinetics = Arrhenius(A=(2.18e+16,'s^-1'), n=0, Ea=(2.9288,'kJ/mol'), T0=(1,'K'), comment="""Estimated using template [R3;Y_rad_out;Ypri_rad_out] for rate rule [R3_SD;Y_rad_out;CdsinglepriDe_rad_out]
Euclidian distance = 2.2360679775
family: Birad_recombination"""),
)

reaction(
    label = 'reaction70',
    reactants = ['[CH2][C]=CC#CC[O](25018)'],
    products = ['C=[C]C=C[C]=C[O](22614)'],
    transitionState = 'TS70',
    kinetics = Arrhenius(A=(2.56742e+09,'s^-1'), n=1.0541, Ea=(193.078,'kJ/mol'), T0=(1,'K'), comment="""Estimated using template [R3Hall;Cd_rad_out_singleDe_Cd;XH_out] for rate rule [R3HJ;Cd_rad_out_singleDe_Cd;XH_out]
Euclidian distance = 1.0
Multiplied by reaction path degeneracy 2.0
family: intra_H_migration"""),
)

reaction(
    label = 'reaction71',
    reactants = ['[CH2]C#CC=[C]C[O](25019)'],
    products = ['C=[C]C=C[C]=C[O](22614)'],
    transitionState = 'TS71',
    kinetics = Arrhenius(A=(2.26683e+09,'s^-1'), n=1.04717, Ea=(150.072,'kJ/mol'), T0=(1,'K'), comment="""Estimated using template [R4Hall;Cd_rad_out_double;XH_out] for rate rule [R4HJ_2;Cd_rad_out_double;XH_out]
Euclidian distance = 1.0
Multiplied by reaction path degeneracy 2.0
family: intra_H_migration"""),
)

reaction(
    label = 'reaction72',
    reactants = ['[CH]=C=CC=[C]C[O](25020)'],
    products = ['C=[C]C=C[C]=C[O](22614)'],
    transitionState = 'TS72',
    kinetics = Arrhenius(A=(3.73886e+06,'s^-1'), n=1.85754, Ea=(151.521,'kJ/mol'), T0=(1,'K'), comment="""Estimated using template [RnH;Cd_rad_out_singleH;XH_out] for rate rule [R6HJ_4;Cd_rad_out_singleH;XH_out]
Euclidian distance = 2.0
Multiplied by reaction path degeneracy 2.0
family: intra_H_migration"""),
)

reaction(
    label = 'reaction73',
    reactants = ['C=[C]C[C][C]=C[O](25021)'],
    products = ['C=[C]C=C[C]=C[O](22614)'],
    transitionState = 'TS73',
    kinetics = Arrhenius(A=(2.94659e+10,'s^-1'), n=0.2847, Ea=(27.8529,'kJ/mol'), T0=(1,'K'), comment="""Estimated using average of templates [Rn;Y_rad;XH_Rrad_De] + [R2radExo;Y_rad;XH_Rrad] for rate rule [R2radExo;Y_rad;XH_Rrad_De]
Euclidian distance = 1.0
Multiplied by reaction path degeneracy 2.0
family: Intra_Disproportionation"""),
)

reaction(
    label = 'reaction74',
    reactants = ['[CH2][CH]C=[C][C]=C[O](25022)'],
    products = ['C=[C]C=C[C]=C[O](22614)'],
    transitionState = 'TS74',
    kinetics = Arrhenius(A=(7.437e+08,'s^-1'), n=1.045, Ea=(63.4002,'kJ/mol'), T0=(1,'K'), comment="""Estimated using an average for rate rule [R3radExo;Y_rad;XH_Rrad]
Euclidian distance = 0
family: Intra_Disproportionation"""),
)

reaction(
    label = 'reaction75',
    reactants = ['C[C]C=[C][C]=C[O](25023)'],
    products = ['C=[C]C=C[C]=C[O](22614)'],
    transitionState = 'TS75',
    kinetics = Arrhenius(A=(1.54267e+09,'s^-1'), n=0.311, Ea=(24.9733,'kJ/mol'), T0=(1,'K'), comment="""Estimated using template [R4;Y_rad;XH_Rrad] for rate rule [R4radEndo;Y_rad;XH_Rrad]
Euclidian distance = 1.0
Multiplied by reaction path degeneracy 3.0
family: Intra_Disproportionation"""),
)

reaction(
    label = 'reaction76',
    reactants = ['C=[C]C=C[C]=C[O](22614)'],
    products = ['C=C=CC=C=C[O](22623)'],
    transitionState = 'TS76',
    kinetics = Arrhenius(A=(1e+10,'s^-1'), n=0, Ea=(0,'kJ/mol'), T0=(1,'K'), comment="""Estimated using an average for rate rule [NOS]
Euclidian distance = 0
family: 1,2-Birad_to_alkene"""),
)

reaction(
    label = 'reaction77',
    reactants = ['C=[C]C=C[C]=C[O](22614)'],
    products = ['C=C=CC1[C]C1[O](25024)'],
    transitionState = 'TS77',
    kinetics = Arrhenius(A=(1e+10,'s^-1'), n=0, Ea=(284.653,'kJ/mol'), T0=(1,'K'), comment="""Estimated using template [R4;multiplebond_intra;radadd_intra_csHCd] for rate rule [R4;carbonylbond_intra_H;radadd_intra_csHCd]
Euclidian distance = 2.0
family: Intra_R_Add_Exocyclic
Ea raised from 281.7 to 284.7 kJ/mol to match endothermicity of reaction."""),
)

network(
    label = '4917',
    isomers = [
        'C=[C]C=C[C]=C[O](22614)',
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
    label = '4917',
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

