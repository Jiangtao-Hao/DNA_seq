function DNA_Seq
    %�õ�����Ƶ��
    %Get_SampleFrequency
    
    %SPSS����
    %����SPSS���ɷַ������õ����ɷ�ϵ��Pcc
    %����SPSS�б𣬵õ��б���ϵ��CDFC�Լ��������Ĵ���ֵ
    
    %�õ�Ԥ��Ƶ��
    %Get_PredictFrequency
    
    %������Ʒ����Ƶ��
    load DNA_PredictFrequency
    load DNA_SampleFrequency
    %���ɷ�ϵ��Pcc
    load Pcc
    %�б���ϵ��
    load CDFC
    %�������Ĵ���ֵ
    load GroupCentroid
    
    %ѡ��DNAƵ��Ϊ��֪���Ǵ���
    Choice = input('ѡ����֪�������м��飬������ 1 \nѡ�������������Ԥ�⣬������ 2');
    if Choice == 1
        DNA_Frequency = DNA_SampleFrequency;
    else
        DNA_Frequency = DNA_PredictFrequency;
    end
    
    
    %������Ʒ���ɷ�Ƶ��
    DNA_MainCom = DNA_Frequency*Pcc;

    %������
    DNA_MainCom(:,8) = 1;
    %�õ�ÿ���������ݵ��б���ֵ
    CDFV = DNA_MainCom * CDFC;
    
    %���������Ʒ������ľ���ֵ����
    for i = 1:size(CDFV,1)
        DF1 = abs(CDFV(i) - GroupCentroid(1));
        DF2 = abs(CDFV(i) - GroupCentroid(2));
        
        if DF1 < DF2
            Result(i) = 1;
        else
            Result(i) = 2;
        end
    end
    
    disp(Result');
end

%�õ�������ƷƵ��
function Get_PredictFrequency
    P1 = 'tttagctcagtccagctagctagtttacaatttcgacaccagtttcgcaccatcttaaatttcgatccgtaccgtaatttagcttagatttggatttaaaggatttagattga';
    P2 = 'tttagtacagtagctcagtccaagaacgatgtttaccgtaacgtqacgtaccgtacgctaccgttaccggattccggaaagccgattaaggaccgatcgaaaggg ';
    P3 = 'cgggcggatttaggccgacggggacccgggattcgggacccgaggaaattcccggattaaggtttagcttcccgggatttagggcccggatggctgggaccc';
    P4 = 'tttagctagctactttagctatttttagtagctagccagcctttaaggctagctttagctagcattgttctttattgggacccaagttcgacttttacgatttagttttgaccgt';
    P5 = 'gaccaaaggtgggctttagggacccgatgctttagtcgcagctggaccagttccccagggtattaggcaaaagctgacgggcaattgcaatttaggcttaggcca';
    P6 = 'gatttactttagcatttttagctgacgttagcaagcattagctttagccaatttcgcatttgccagtttcgcagctcagttttaacgcgggatctttagcttcaagctttttac ';
    P7 = 'ggattcggatttacccggggattggcggaacgggacctttaggtcgggacccattaggagtaaatgccaaaggacgctggtttagccagtccgttaaggcttag';
    P8 = 'tccttagatttcagttactatatttgacttacagtctttgagatttcccttacgattttgacttaaaatttagacgttagggcttatcagttatggattaatttagcttattttcga';
    P9 = 'ggccaattccggtaggaaggtgatggcccgggggttcccgggaggatttaggctgacgggccggccatttcggtttagggagggccgggacgcgttagggc';
    P10 ='cgctaagcagctcaagctcagtcagtcacgtttgccaagtcagtaatttgccaaagttaaccgttagctgacgctgaacgctaaacagtattagctgatgactcgta';
    P11 = 'ttaaggacttaggctttagcagttactttagtttagttccaagctacgtttacgggaccagatgctagctagcaatttattatccgtattaggcttaccgtaggtttagcgt';
    P12 = 'gctaccgggcagtctttaacgtagctaccgtttagtttgggcccagccttgcggtgtttcggattaaattcgttgtcagtcgctctrtgggtttagtcattcccaaaagg';
    P13 = 'cagttagctgaatcgtttagccatttgacgtaaacatgattttacgtacgtaaattttagccctgacgtttagctaggaatttatgctgacgtagcgatcgactttagcac';
    P14 = 'cggttagggcaaaggttggatttcgacccagggggaaagcccgggacccgaacccagggctttagcgtaggctgacgctaggcttaggttggaacccggaaa';
    P15 = 'gcggaagggcgtaggtttgggatgcttagccgtaggctagctttcgacacgatcgattcgcaccacaggataaaagttaagggaccggtaagtcgcggtagcc';
    P16 = 'ctagctacgaacgctttaggcgcccccgggagtagtcgttaccgttagtatagcagtcgcagtcgcaattcgcaaaagtccccagctttagccccagagtcgacg';
    P17 = 'gggatgctgacgctggttagctttaggcttagcgtagctttagggccccagtctgcaggaaatgcccaaaggaggcccaccgggtagatgccasagtgcaccgt';
    P18 = 'aacttttagggcatttccagttttacgggttattttcccagttaaactttgcaccattttacgtgttacgatttacgtataatttgaccttattttggacactttagtttgggttac';
    P19 = 'ttagggccaagtcccgaggcaaggaattctgatccaagtccaatcacgtacagtccaagtcaccgtttgcagctaccgtttaccgtacgttgcaagtcaaatccat';
    P20 = 'ccattagggtttatttacctgtttattttttcccgagaccttaggtttaccgtactttttaacggtttacctttgaaatttttggactagcttaccctggatttaacggccagttt';
    
    ATCG_AT = ones(20,5);
    AAF = ones(20,21);
    for i = 1:20
        DNA_name = ['P',num2str(i)];
        
        ATCG_AT(i,:) = BaseFrequency(eval(DNA_name));
        AAF(i,:) = AminoAcidFrequency_Best(eval(DNA_name));
    end

    DNA_PredictFrequency = [ATCG_AT,AAF];
    
    save DNA_PredictFrequency
end

%�õ���֪��ƷƵ��
function Get_SampleFrequency
A1 = 'aggcacggaaaaacgggaataacggaggaggacttggcacggcattacacggaggacgaggtaaaggaggcttgtctacggccggaagtgaagggggatatgaccgcttgg';
A2 = 'cggaggacaaacgggatggcggtattggaggtggcggactgttcggggaattattcggtttaaacgggacaaggaaggcggctggaacaaccggacggtggcagcaaagga';
A3 = 'gggacggatacggattctggccacggacggaaaggaggacacggcggacatacacggcggcaacggacggaacggaggaaggagggcggcaatcggtacggaggcggcgga';
A4 = 'atggataacggaaacaaaccagacaaacttcggtagaaatacagaagcttagatgcatatgttttttaaataaaatttgtattattatggtatcataaaaaaaggttgcga';
A5 = 'cggctggcggacaacggactggcggattccaaaaacggaggaggcggacggaggctacaccaccgtttcggcggaaaggcggagggctggcaggaggctcattacggggag';
A6 = 'atggaaaattttcggaaaggcggcaggcaggaggcaaaggcggaaaggaaggaaacggcggatatttcggaagtggatattaggagggcggaataaaggaacggcggcaca';
A7 = 'atgggattattgaatggcggaggaagatccggaataaaatatggcggaaagaacttgttttcggaaatggaaaaaggactaggaatcggcggcaggaaggatatggaggcg';
A8 = 'atggccgatcggcttaggctggaaggaacaaataggcggaattaaggaaggcgttctcgcttttcgacaaggaggcggaccataggaggcggattaggaacggttatgagg';
A9 = 'atggcggaaaaaggaaatgtttggcatcggcgggctccggcaactggaggttcggccatggaggcgaaaatcgtgggcggcggcagcgctggccggagtttgaggagcgcg';
A10 ='tggccgcggaggggcccgtcgggcgcggatttctacaagggcttcctgttaaggaggtggcatccaggcgtcgcacgctcggcgcggcaggaggcacgcgggaaaaaacg';

B1 = 'gttagatttaacgttttttatggaatttatggaattataaatttaaaaatttatattttttaggtaagtaatccaacgtttttattactttttaaaattaaatatttatt';
B2 = 'gtttaattactttatcatttaatttaggttttaattttaaatttaatttaggtaagatgaatttggttttttttaaggtagttatttaattatcgttaaggaaagttaaa';
B3 = 'gtattacaggcagaccttatttaggttattattattatttggattttttttttttttttttttaagttaaccgaattattttctttaaagacgttacttaatgtcaatgc';
B4 = 'gttagtcttttttagattaaattattagattatgcagtttttttacataagaaaatttttttttcggagttcatattctaatctgtctttattaaatcttagagatatta';
B5 = 'gtattatatttttttatttttattattttagaatataatttgaggtatgtgtttaaaaaaaatttttttttttttttttttttttttttttttaaaatttataaatttaa';
B6 = 'gttatttttaaatttaattttaattttaaaatacaaaatttttactttctaaaattggtctctggatcgataatgtaaacttattgaatctatagaattacattattgat';
B7 = 'gtatgtctatttcacggaagaatgcaccactatatgatttgaaattatctatggctaaaaaccctcagtaaaatcaatccctaaacccttaaaaaacggcggcctatccc';
B8 = 'gttaattatttattccttacgggcaattaattatttattacggttttatttacaattttttttttttgtcctatagagaaattacttacaaaacgttattttacatactt';
B9 = 'gttacattatttattattatccgttatcgataattttttacctcttttttcgctgagtttttattcttactttttttcttctttatataggatctcatttaatatcttaa';
B10 = 'gtatttaactctctttactttttttttcactctctacattttcatcttctaaaactgtttgatttaaacttttgtttctttaaggattttttttacttatcctctgttat';

ATCG_AT = ones(20,5);
AAF = ones(20,21);
for i = 1:20
    if i <= 10
        DNA_name = ['A',num2str(i)];
    else
        DNA_name = ['B',num2str(i-10)];
    end
    
    ATCG_AT(i,:) = BaseFrequency(eval(DNA_name));
    AAF(i,:) = AminoAcidFrequency_Best(eval(DNA_name));
end

DNA_SampleFrequency = [ATCG_AT,AAF];
save DNA_SampleFrequency

%��DNA_SampleFrequency��ΪExcel�ļ�������SPSS����
[m,n] = size(Data);
DNA_SampleFrequency_cell = mat2cell(DNA_SampleFrequency,ones(m,1),ones(n,1));

%ȡ��
for i = 1:21
    Title_AA(i) = i;
end
Title_AA = num2str(Title_AA);
Title_AA = split(Title_AA);

Title_ATCG = {'A','T','C','G','A+T'};
Title = [Title_ATCG,Title_AA'];


Result = [Title;DNA_SampleFrequency_cell];

xlswrite('DNA_SampleFrequency.xls',Result);
end

%�õ����Ƶ��
function ATCG_AT = BaseFrequency(DNA)
ATCG_AT = zeros(1,5);

for i=1:length(DNA)
    if DNA(i) == 'a'
        ATCG_AT(1) = ATCG_AT(1) + 1;
    end
    
    if DNA(i) == 't'
        ATCG_AT(2) = ATCG_AT(2) + 1;
    end
    
    if DNA(i) == 'c'
        ATCG_AT(3) = ATCG_AT(3) + 1;
    end
    
    if DNA(i) == 'g'
        ATCG_AT(4) = ATCG_AT(4) + 1;
    end
end
ATCG_AT(5) = ATCG_AT(1) + ATCG_AT(2);

ATCG_AT = ATCG_AT./length(DNA);

end

%�õ�������Ƶ��
function AAF = AminoAcidFrequency(DNA)
AAF = zeros(1,21);%21�ְ�����

Condon_1 = {'aaa','aag'};
Condon_2 = {'aat','aac','gaa','gag','gat','gac'};
Condon_3 = {'aga','agg','agt','agc','tca','tcg'};
Condon_4 = {'ata','atg'};
Condon_5 = {'aca','acg'};
Condon_6 = {'acc'};
Condon_7 = {'gga','ggg','ggt','ggc'};
Condon_8 = {'gat','gtg'};
Condon_9 = {'gtt','gtc'};
Condon_10 = {'gca','gcg','gct','gcc','tct','tcc'};
Condon_11 = {'taa','tag','tat'};
Condon_12 = {'tga','tgg','tgt','tgc'};
Condon_13 = {'tta','ttg'};
Condon_14 = {'ttt','ttc'};
Condon_15 = {'caa','cag','cat','cac'};
Condon_16 = {'cga','cgg','cgt','cgc'};
Condon_17 = {'cta','ctg'};
Condon_18 = {'ctt','ctc'};
Condon_19 = {'cca','ccg','cct','ccc'};
Condon_20 = {'tac'};
Condon_21 = {'att','atc','act'};


for i = 1:length(DNA)/3
    Condon = DNA(3*i-2:3*i);
    
    for j = 1:21
        AAF_Name = ['Condon_',num2str(j)];
        
        MatchInformation = strcmp(Condon,eval(AAF_Name));
        
        if any(MatchInformation)
            AAF(j) = AAF(j) + 1;
        end
    end
end

AAF = AAF./(length(DNA)/3);
end

%ѡ�����Žض���ʽ
function AAF_Best = AminoAcidFrequency_Best(DNA)
    P = [2 6 6 2 2 1 4 2 2 6 3 4 2 2 4 4 2 2 4 1 3];
    P = P./64;
    
    AAF_1 = AminoAcidFrequency(DNA);%��һ�ֽض���ʽ
    AAF_2 = AminoAcidFrequency(DNA(2:end));%�ڶ��ֽض���ʽ
    AAF_3 = AminoAcidFrequency(DNA(3:end));%�����ֽض���ʽ
    
    E = P*[AAF_1' AAF_2' AAF_3'];
    [MaxE,MaxOrder] = max(E);
    
    AAF_Best = eval(['AAF_',num2str(MaxOrder)]);
end
