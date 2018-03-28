%% PIFFER METHOD OF ANALYZING RACIAL DISTIBUTION OF SNPS CORRELATED WITH INTELLIGENCE

% Creating intial files / folders / structures / tables

% Create two folders in MATLAB folder. One contains the allele info (rs.INFO files), and other the
% distiribution in population (rs.FORMAT files). Rename rs.format files as SNPnumber.txt (i.e SNP1.txt)
% and rs.info as rs.INFO.txt (i.e rs236330.INFO.txt). It does not matter how you number the SNPnumber files but start
% from 1 upwards. Add both folders to path. The samples file rename as samples.txt. Create an excel file that in one
% coulumn has SNP name (i.e rs236330), second column has effect allele (i.e A) column 3 is beta, column 4 is
% Z statistic, and column 5 is effect direction. For each SNP only one value needs to be filled. Can only work for
% SNPs with 2 distinct alleles. For 3 alleles + will write code later as not common issue.

tic; SNP_projectdata=struct(); popref=struct();
IQmat={'IQ',83,85,81,91,99,105,105,83.5,71,101,100,82,62,97,82,105,99.4,74,64,88,85,84,83.5,79,99,71};
[raw1,raw2,SNP_efdata] = xlsread('EffectData'); % excel import for effectdata here
SNP_format = dir('SNP_format1'); Allele_info = dir('allele_info1'); % rename folders here
if strcmp(SNP_format(3).name, '.DS_Store')
    SNP_format(3) = [];
end
if strcmp(Allele_info(3).name, '.DS_Store')
    Allele_info(3) = [];
end
SNP_number = size(SNP_format); SNP_number = (SNP_number (1,1))-3; A = importdata ('samples.txt');
T = [{'SNP','ANC','DER'},A.textdata]; SNP_projectdata.rawdata.headers = T;
SNP_projectdata.calcdata_subpop.headers = A.textdata(3:29); SNP_POStableFA = A.textdata(4:29);
SNP_projectdata.calcdata_superpop.headers = {'global','AFR','EAS','EUR','SAS','AMR'};
SNP_rawtable = SNP_projectdata.rawdata.headers; SNP_POStable = {'AFR','EAS','EUR','SAS','AMR'};
SNP_rawtable(SNP_number+4,6:32) = IQmat; disp ' '; 
disp '1 - SNP information loaded'; disp(cat(2,'Number of SNPs entered = ',num2str(SNP_number))); toc; disp ' ';

%% Totalling populations

for run1=1:27
    SNP_rawtable(2,(run1+5)) = {A.data(run1+2)}; SNP_projectdata.rawdata.poptotal(run1) = A.data(run1+2);
    SNP_projectdata.calcdata_subpop.subtotal(run1) = A.data(run1+2);
end

popref.p1 = 1; popref.p2 = [2,3,10,14,19,20,27]; popref.p3 = [7,8,17,18,5]; popref.p4 = [6,11,12,15,26];
popref.p5 = [4,13,16,23,25]; popref.p6 = [9,21,22,24]; popref.pom1 = 1; popref.pom2 = [2,3,9,12,16,17,24];
popref.pom3 = [6,7,14,15]; popref.pom4 = [5,10,11,13,23]; popref.pom5 = [4,20,22]; popref.pom6 = [8,18,19,21];

for run2 = 1:6
    superef=size(popref.(cat(2,'p',(num2str(run2))))); superef=superef(1,2);
    
    for run3 = 1:superef
        poprefnum=(popref.(cat(2,'p',(num2str(run2))))); poprefnum=poprefnum(1,run3);
        if run3 == 1
            SNP_projectdata.calcdata_superpop.superpoptotal(run2) = SNP_projectdata.rawdata.poptotal(poprefnum);
        else
            SNP_projectdata.calcdata_superpop.superpoptotal(run2) = ...
                SNP_projectdata.calcdata_superpop.superpoptotal(run2) + SNP_projectdata.rawdata.poptotal(poprefnum);
        end
    end
    
end

disp '2 - Populations totalled'; disp (cat(2,'Total = ',mat2str(cell2mat(SNP_rawtable(2,6)))));
disp(cat(2,cat(2,'AFR = ',mat2str(SNP_projectdata.calcdata_superpop.superpoptotal(2))),' (ACB,ASW,ESN,GWD,LWK,MSL,YRI)'));
disp(cat(2,cat(2,'EAS = ',mat2str(SNP_projectdata.calcdata_superpop.superpoptotal(3))),' (CDX,CHB,CHS,JPT,KHV)'));
disp(cat(2,cat(2,'EUR = ',mat2str(SNP_projectdata.calcdata_superpop.superpoptotal(4))),' (CEU,FIN,GBR,IBS,TSI)'));
disp(cat(2,cat(2,'SAS = ',mat2str(SNP_projectdata.calcdata_superpop.superpoptotal(5))),' (BEB,GIH,ITU,PJL,STU)'));
disp(cat(2,cat(2,'AMR = ',mat2str(SNP_projectdata.calcdata_superpop.superpoptotal(6))),' (CLM,MXL,PEL,PUR)')); toc; disp ' ';

%% Input / format / manipulate raw SNP data (subpop)

for run4 = 1:2
    
    for run5=1:SNP_number
        A = importdata(cat(2,(cat(2,'SNP',(num2str(run5)))),'.txt'));
        SNP_projectdata.rawdata.(cat(2,'SNP',(num2str(run5)))) = A.data;
        
        for run6=1:27
            SNP_projectdata.calcdata_subpop.(cat(2,'SNP',(num2str(run5)))).ANC_allele_count(run6) = ...
                A.data(run6+2);
            SNP_projectdata.calcdata_subpop.(cat(2,'SNP',(num2str(run5)))).ALT_allele_count(run6) = ...
                cell2mat(SNP_rawtable(2,run6+5)) - A.data(run6+2);
            SNP_projectdata.calcdata_subpop.(cat(2,'SNP',(num2str(run5)))).ANC_allele_freq(run6) = ...
                (A.data(run6+2))*100 / cell2mat(SNP_rawtable(2,run6+5));
            SNP_projectdata.calcdata_subpop.(cat(2,'SNP',(num2str(run5)))).ALT_allele_freq(run6) = ...
                (cell2mat(SNP_rawtable(2,run6+5)) - A.data(run6+2))*100 / cell2mat(SNP_rawtable(2,run6+5));
        end
        
        if run4 == 2
            charlim = size(cell2mat(SNP_rawtable((run5+2),3))); charlim = charlim(1,2);
            if charlim == 9 || charlim > 9
                
                for run7=1:27
                    SNP_projectdata.calcdata_subpop.(cat(2,'SNP',(num2str(run5)))).ALT_allele_count(run7) = ...
                        A.data(run7+2);
                    SNP_projectdata.calcdata_subpop.(cat(2,'SNP',(num2str(run5)))).ANC_allele_count(run7) = ...
                        cell2mat(SNP_rawtable(2,run7+5)) - A.data(run7+2);
                    SNP_projectdata.calcdata_subpop.(cat(2,'SNP',(num2str(run5)))).ALT_allele_freq(run7) = ...
                        (A.data(run7+2))*100 / cell2mat(SNP_rawtable(2,run7+5));
                    SNP_projectdata.calcdata_subpop.(cat(2,'SNP',(num2str(run5)))).ANC_allele_freq(run7) = ...
                        (cell2mat(SNP_rawtable(2,run7+5)) - A.data(run7+2))*100 / cell2mat(SNP_rawtable(2,run7+5));
                end
                
            end
        end
        
        for run8=1:29
            SNP_rawtable((run5+2),(run8+3)) = {A.data(run8)};
        end
        
    end
    
    if run4 == 2
        disp '4 - Allele count / frequencies calculated for sub populations'; toc; disp ' ';
        break
    end
    
    %% Allele info prompts
    
    disp '3 - Answer following prompts with (y/n)'; SNP_projectdata.rawdata.Allele_info = Allele_info;
    prompt = 'Is the impact of alleles on trait known? (enter y/n): '; Known = input(prompt,'s');
    if Known == 'y'
        ValueKnown = 1; prompt = 'Use excel import(y) or manual input(n)? (enter y/n): '; prevdir = input(prompt,'s');
        if prevdir == 'y'
            ValueKnown1 = 'y';
        else
            ValueKnown1 = 'n';
        end
    else
        ValueKnown = 0;
    end
    prompt = 'Show bar graphs which will take a lot longer? (enter y/n): '; BARon = input(prompt,'s'); toc; disp ' ';
    
    for run9=1:SNP_number
        Allele = Allele_info(run9+2).name; T = readtable (Allele); AA = cell2mat(table2cell(T(1,5)));
        if upper(AA(1:end-3)) > 0
            AA = upper(AA(1:end-3));
        end
        if AA == upper(cell2mat(table2cell(T(1,3))))
            aANC = table2cell(T(1,3)); aDER = table2cell(T(1,4));
        else
            aANC = table2cell(T(1,4)); aDER = table2cell(T(1,3));
        end
        REF = cell2mat(table2cell(T(1,4))); REF = upper(REF);
        
        for run10=1:SNP_number
            if cell2mat(table2cell(T(1,2))) == cell2mat(SNP_rawtable(run10+2,5))
                SNP_rawtable((run10+2),1) = {Allele(1:end-9)}; n=run10+2;
            end
        end
        
        for run11=1:(size(SNP_efdata,1)-1)
            rsc = cell2mat(SNP_rawtable(n,1)); rsc1 = cell2mat(SNP_efdata(run11+1,1));
            if size(rsc,2) == size(rsc1,2)
                if rsc == rsc1
                    Aref = (cell2mat(SNP_efdata(run11+1,2)));
                    if double(isnan(cell2mat(SNP_efdata(run11+1,5)))) == 0
                        Adir = (cell2mat(SNP_efdata(run11+1,5))); impact = Adir;
                    else
                        if double(isnan(cell2mat(SNP_efdata(run11+1,3)))) == 0
                            Ab = cell2mat(SNP_efdata(run11+1,3));
                            if Ab < 0
                                impact = '-';
                            else
                                impact = '+';
                            end
                        else
                            if double(isnan(cell2mat(SNP_efdata(run11+1,4)))) == 0
                                Az = cell2mat(SNP_efdata(run11+1,4));
                                if Az < 0
                                    impact = '-';
                                else
                                    impact = '+';
                                end
                            end
                        end
                    end
                end
            end
        end
        
        if  REF == cell2mat(aANC)
            aANC = {(cat(2,cell2mat(aANC),' (ref)'))};
            if ValueKnown == 1
                if ValueKnown1 == 'y'
                    if Aref == REF
                        aANC = {cat(2,(cat(2,cell2mat(aANC),' ')),impact)};
                        if  impact == '+'
                            aDER = {cat(2,(cat(2,cell2mat(aDER),' ')),'-')};
                        else
                            aDER = {cat(2,(cat(2,cell2mat(aDER),' ')),'+')};
                        end
                    else
                        aDER = {cat(2,(cat(2,cell2mat(aDER),' ')),impact)};
                        if  impact == '+'
                            aANC = {cat(2,(cat(2,cell2mat(aANC),' ')),'-')};
                        else
                            aANC = {cat(2,(cat(2,cell2mat(aANC),' ')),'+')};
                        end
                    end
                end
                if ValueKnown1 == 'n'
                    prompt = cat(2,(cat(2,(cat(2,'Is allele ',REF)),cat(2,' for SNP ',Allele(1:end-9)))),...
                        ' positive or deleterious? (enter: +/-): '); impact = input(prompt,'s');
                    aANC = {cat(2,(cat(2,cell2mat(aANC),' ')),impact)};
                    if  impact == '+'
                        aDER = {cat(2,(cat(2,cell2mat(aDER),' ')),'-')};
                    else
                        aDER = {cat(2,(cat(2,cell2mat(aDER),' ')),'+')};
                    end
                end
            else
                impact = rand(1);
                if impact > 0.5
                    aDER = {cat(2,(cat(2,cell2mat(aDER),' ')),'+')}; aANC = {cat(2,(cat(2,cell2mat(aANC),' ')),'-')};
                else
                    aDER = {cat(2,(cat(2,cell2mat(aDER),' ')),'-')}; aANC = {cat(2,(cat(2,cell2mat(aANC),' ')),'+')};
                end
            end
        else
            aDER = {(cat(2,cell2mat(aDER),' (ref)' ))};
            if ValueKnown == 1
                if ValueKnown1 == 'y'
                    if Aref == REF
                        aDER = {cat(2,(cat(2,cell2mat(aDER),' ')),impact)};
                        if  impact == '+'
                            aANC = {cat(2,(cat(2,cell2mat(aANC),' ')),'-')};
                        else
                            aANC = {cat(2,(cat(2,cell2mat(aANC),' ')),'+')};
                        end
                    else
                        aANC = {cat(2,(cat(2,cell2mat(aANC),' ')),impact)};
                        if  impact == '+'
                            aDER = {cat(2,(cat(2,cell2mat(aDER),' ')),'-')};
                        else
                            aDER = {cat(2,(cat(2,cell2mat(aDER),' ')),'+')};
                        end
                    end
                end
                if ValueKnown1 == 'n'
                    prompt = cat(2,(cat(2,(cat(2,'is allele ',REF)),cat(2,' for SNP ',Allele(1:end-9)))),...
                        ' positive or deleterious (enter: +/-): '); impact = input(prompt,'s');
                    aDER = {cat(2,(cat(2,cell2mat(aDER),' ')),impact)};
                    if  impact == '+'
                        aANC = {cat(2,(cat(2,cell2mat(aANC),' ')),'-')};
                    else
                        aANC = {cat(2,(cat(2,cell2mat(aANC),' ')),'+')};
                    end
                end
            else
                impact = rand(1);
                if impact > 0.5
                    aDER = {cat(2,(cat(2,cell2mat(aDER),' ')),'+')}; aANC = {cat(2,(cat(2,cell2mat(aANC),' ')),'-')};
                else
                    aDER = {cat(2,(cat(2,cell2mat(aDER),' ')),'-')}; aANC = {cat(2,(cat(2,cell2mat(aANC),' ')),'+')};
                end
            end
        end
        
        for run12=1:SNP_number
            if cell2mat(table2cell(T(1,2))) == cell2mat(SNP_rawtable(run12+2,5))
                SNP_rawtable((run12+2),2) = aANC; SNP_rawtable((run12+2),3) = aDER;
                popref.(cat(2,'r',(num2str(run12)))) = REF;
            end
        end
        
    end
    
end

%% Input / format / manipulate raw SNP data (superpop)

for run13=1:SNP_number
    
    for run14=1:4
        stat1 = {'ALT_allele_count','ANC_allele_count','ALT_allele_freq','ANC_allele_freq'}; stat2 = stat1(run14);
        
        for run15 = 1:6
            superef=size(popref.(cat(2,'p',(num2str(run15))))); superef=superef(1,2);
            if run14 == 1 || run14 == 2 || run14 == 3 || run14 == 4
                
                for run16 = 1:superef
                    poprefnum=(popref.(cat(2,'p',(num2str(run15))))); poprefnum=poprefnum(1,run16);
                    if run16 == 1
                        SNP_projectdata.calcdata_superpop.(cat(2,'SNP',(num2str(run13)))).(char(stat2))(run15) = ...
                            SNP_projectdata.calcdata_subpop.(cat(2,'SNP',(num2str(run13)))).(char(stat2))(poprefnum);
                        statgroup = ...
                            SNP_projectdata.calcdata_subpop.(cat(2,'SNP',(num2str(run13)))).(char(stat2))(poprefnum);
                    else
                        SNP_projectdata.calcdata_superpop.(cat(2,'SNP',(num2str(run13)))).(char(stat2))(run15) = ...
                            SNP_projectdata.calcdata_superpop.(cat(2,'SNP',(num2str(run13)))).(char(stat2))(run15) + ...
                            SNP_projectdata.calcdata_subpop.(cat(2,'SNP',(num2str(run13)))).(char(stat2))(poprefnum);
                        statgroup = [statgroup,...
                            SNP_projectdata.calcdata_subpop.(cat(2,'SNP',(num2str(run13)))).(char(stat2))(poprefnum)];
                    end
                end
                
            end
            if run14 == 3 || run14 == 4
                stat3 = stat1(run14-2);
                SNP_projectdata.calcdata_superpop.(cat(2,'SNP',(num2str(run13)))).(char(stat2))(run15) = ...
                    SNP_projectdata.calcdata_superpop.(cat(2,'SNP',(num2str(run13)))).(char(stat3))(run15)*100 /...
                    SNP_projectdata.calcdata_superpop.superpoptotal(run15);
            else
            end
            dev=std(statgroup);
            SNP_projectdata.calcdata_superpop.(cat(2,'SNP',(num2str(run13)))).(cat(2,(char(stat2)),'dev'))(run15) = dev;
        end
        
    end
end

disp '5 - Allele count / frequencies / stdev calculated for super populations'; toc; disp ' ';

%% Raw SNP data table / REF allele count table

figure('name','RawSNPdata')
set(gcf, 'Position', [0 1000 1450 270])
Tableout1=uitable('Position',[0 0 1450 270],'Data',SNP_rawtable); disp '6 - Raw table generated'; toc; disp ' ';

%% Positive allele frequencies in populations by super population bar graphs

if mod(SNP_number,9) == 0
    pl = floor(SNP_number/9);
else
    pl = ceil(SNP_number/9);
end
if pl<1
    pl=1;
else
end
if SNP_number > 27
    pl=3;
end
n = 0; m = 1;
if BARon == 'y'
    figure('name','percent positive allele in populations by super population')
    set(gcf, 'Position', [1000 105 1350 680])
    X = categorical({'global','African','East Asian','European','South Asian','Americas'});
end

for run17 = 1:SNP_number
    if BARon == 'y'
        if SNP_number > 36
            if n > 35
                figure('name','percent positive allele in populations by super population')
                set(gcf, 'Position', [1000 105 1350 680])
                n=0; m=m+1;
            end
        end
        n=n+1; ax = subplot(pl+1,9,n);
    end
    ef = SNP_rawtable((run17+2),2); ef = cell2mat(ef); ef = ef(end);
    if ef == '+'
        Ypos = SNP_projectdata.calcdata_superpop.(cat(2,'SNP',(num2str(run17)))).ANC_allele_freq;
        YposFA = SNP_projectdata.calcdata_subpop.(cat(2,'SNP',(num2str(run17)))).ANC_allele_freq;
        if BARon == 'y'
            titleRS = strjoin(cat(2,(SNP_rawtable(run17+2)),SNP_rawtable((run17+2),2)));
            err = SNP_projectdata.calcdata_superpop.(cat(2,'SNP',(num2str(run17)))).ANC_allele_freqdev;
        end
    else
        Ypos = SNP_projectdata.calcdata_superpop.(cat(2,'SNP',(num2str(run17)))).ALT_allele_freq;
        YposFA = SNP_projectdata.calcdata_subpop.(cat(2,'SNP',(num2str(run17)))).ALT_allele_freq;
        if BARon == 'y'
            titleRS = strjoin(cat(2,(SNP_rawtable(run17+2)),SNP_rawtable((run17+2),3)));
            err = SNP_projectdata.calcdata_superpop.(cat(2,'SNP',(num2str(run17)))).ALT_allele_freqdev;
        end
    end
    SNP_POStable(run17+1,:) = num2cell(Ypos(2:6)); SNP_POStableFA(run17+1,:) = num2cell(YposFA(2:27));
    if BARon == 'y'
        hBar = bar(ax,X,Ypos); title(titleRS); ylim([0 100]);
        hold on
        hBar.FaceColor = 'flat'; hBar.CData(6,:) = [.5 0 .5]; errorbar(ax,X,Ypos,err,'.');
    end
end

disp '7 - Bargraph of frequency positive allele by super population'; 
if BARon == 'y'
disp (cat(2,'Number of pages with bar graphs (1 page contains max 36 bar graphs) = ',num2str(m)))
else
    disp 'Bar graphs not displayed';
end
toc; disp ' ';

%% Anova boxplot for positive allele frequencies w/ p stat by super population

if SNP_number > 27 || BARon == 'n'
    figure('name','Anova boxplot')
    set(gcf, 'Position', [1000 480 750 400])
end
hold on
ax = subplot(pl+1,9,[(pl*9)+5,(pl*9)+9]);
if SNP_number > 27 || BARon == 'n'
    ax = subplot(1,1,1);
end

Anovadata = [(cell2mat(SNP_POStable(2:SNP_number+1,1)))',(cell2mat(SNP_POStable(2:SNP_number+1,2)))',...
    (cell2mat(SNP_POStable(2:SNP_number+1,3)))',(cell2mat(SNP_POStable(2:SNP_number+1,4)))',...
    (cell2mat(SNP_POStable(2:SNP_number+1,5)))'];
Anovacat = [repmat({'AFR'},[1 SNP_number]),repmat({'EAS'},[1 SNP_number]),repmat({'EUR'},[1 SNP_number]),...
    repmat({'SAS'},[1 SNP_number]),repmat({'AMR'},[1 SNP_number])];
[P_anova,Ptbl,Pstat] = anova1(Anovadata,Anovacat,'off'); dim = [0.91 0 0.3 0.3];
boxplot(ax,Anovadata,Anovacat); ylim([0 100]); title('Anova boxplot for positive allele');
annotation('textbox',dim,'String',cat(2,'P = ',num2str(P_anova)),'FitBoxToText','on');
disp '8 - Anova boxplot of frequency positive allele'; disp (cat(2,'P_anova = ',num2str(P_anova))); toc; disp ' ';

%% FA positive allele / polygenic score and regression

SNP_POStableFAc1=cell2mat(SNP_POStableFA(2:SNP_number+1,:)'); SNP_POStableFAc2 = SNP_POStableFAc1;
a21 = cov(SNP_POStableFAc2); a22 = eig(a21); a23 = SNP_number; a24 = 0; a25 = 0;

for run18 = 1:SNP_number
    a24 = a24 +1;
    if a22(a24) < 0.0000001
        SNP_POStableFAc2(:,1)=[]; a21 = cov(SNP_POStableFAc2); a22 = eig(a21); a23 = a23-1; a24 = 0;
        if run18 == 1
            a25 = 1;
        end
    else
        break
    end
end

[Fa_loading,Fa_psi,Fa_T,Fa_stats,Fa_score]=factoran(SNP_POStableFAc2,1,'scores','regression');
PolygenicScore = (SNP_rawtable(1,7:32))';

for run19 = 1:26
    PolygenicScore(run19,2)={sum(SNP_POStableFAc1(run19,:))/SNP_number};
    PolygenicScore(run19,4)=IQmat(run19+1); PolygenicScore(run19,5)={Fa_score(run19)};
end

Fa_score=PolygenicScore(:,[1 5]); PolygenicScoreOM=PolygenicScore; PolygenicScoreOM([4 12 15],:)=[];
X2=cell2mat(PolygenicScore(:,2)); Y2=cell2mat(PolygenicScore(:,4));
pcoef = polyfit(X2,Y2,1); pfit = polyval(pcoef,X2); presid = Y2 - pfit; SSresid = sum(presid.^2);
SStotal = (length(Y2)-1) * var(Y2); rsq = 1 - SSresid/SStotal; R_polygenicscore = rsq^0.5;
X3=cell2mat(PolygenicScoreOM(:,2)); Y3=cell2mat(PolygenicScoreOM(:,4));
pcoef = polyfit(X3,Y3,1); pfit1 = polyval(pcoef,X3); presid = Y3 - pfit1; SSresid = sum(presid.^2);
SStotal = (length(Y3)-1) * var(Y3); rsq = 1 - SSresid/SStotal; R_polygenicscoreOM = rsq^0.5;
disp '9 - Factor score / polygenic score calculated by sub population'; 
if a25 == 1
disp 'Correlation Matix produced negative eig valaues'; disp 'Removed entries to resolve matrix';
end
toc; disp ' ';

%% Positive allele frequencies in populations by sub population bar graphs

col=struct(); n=0; m=1;
col.c1=[.5 0 .5]; col.c2=[0 .5 .5]; col.c3=[.5 .5 0]; col.c4=[1 0 0]; col.c5=[0 0 1];col.c6=[0 1 0];
if mod(SNP_number,4) == 0
    pl1 = floor(SNP_number/4); pl = pl1+3;
else
    pl1 = ceil(SNP_number/4); pl = pl1+3;
end
if pl<1
    pl=2;
end
if SNP_number > 28
    pl = 8;
    if SNP_number > 32
        pl = 9;
    end
end
if BARon == 'n'
else
    figure('name','percent positive allele in populations by sub population')
    set(gcf,'Position',[0 0 1950 750])
    X = categorical(SNP_projectdata.calcdata_subpop.headers);
    
    for run20 = 1:SNP_number
        if SNP_number > 36
            if n > 35
                figure('name','percent positive allele in populations by sub population')
                set(gcf,'Position',[0 0 1950 750])
                n=0; m=m+1;
            end
        end
        n=n+1; ax = subplot(pl,4,n);
        ef = SNP_rawtable((run20+2),2); ef = cell2mat(ef); ef = ef(end);
        if ef == '+'
            Y = SNP_projectdata.calcdata_subpop.(cat(2,'SNP',(num2str(run20)))).ANC_allele_freq;
            titleRS = strjoin(cat(2,(SNP_rawtable(run20+2)),SNP_rawtable((run20+2),2)));
        else
            Y = SNP_projectdata.calcdata_subpop.(cat(2,'SNP',(num2str(run20)))).ALT_allele_freq;
            titleRS = strjoin(cat(2,(SNP_rawtable(run20+2)),SNP_rawtable((run20+2),3)));
        end
        hold on
        
        for run21 = 1:6
            superef=size(popref.(cat(2,'p',(num2str(run21))))); superef=superef(1,2);
            color = col.(cat(2,'c',(num2str(run21)))); poprefnum=(popref.(cat(2,'p',(num2str(run21)))));
            bar(ax,X(poprefnum),Y(poprefnum),'FaceColor',color); title(titleRS);ylim([0 100]); yticks([0 50 100]);
        end
        
    end
    
end
disp '10 - Bargraph of frequency positive allele by sub population';
if BARon == 'y'
disp (cat(2,'Number of pages with bar graphs (1 page contains max 36 bar graphs) = ',num2str(m)))
else
    disp 'Bar graphs not displayed';
end
toc; disp ' ';

%% Scatter plot for polygenic score w/ regression by sub populations

if SNP_number > 28 || BARon == 'n'
    figure('name','Regression analysis of National IQ vs Polygenic Score')
    set(gcf, 'Position', [0 0 1450 300])
end
if SNP_number > 28 || BARon == 'n'
    subplot(1,2,1);
else
    plotnum = [((pl1+1)*4)+1,((pl1+1)*4)+2,((pl1+2)*4)+1,((pl1+2)*4)+2];
    subplot(pl,4,plotnum);
end
plotname = {'AFR','EAS','EUR','SAS','AMR'};
hold on

for run22 = 2:6
    poprefnum=(popref.(cat(2,'p',(num2str(run22)))))-1; color = col.(cat(2,'c',(num2str(run22))));
    plot(X2(poprefnum),Y2(poprefnum),'x','MarkerSize',12,'MarkerFaceColor',color,'MarkerEdgeColor',color,...
        'DisplayName',cell2mat(plotname(run22-1)));
end

ylim([60 110]); title({'';'Average frequency of positive allele in population vs IQ'})
annotation('textbox',[0.03 0 0.3 0.3],'String',cat(2,'r_1 (all) = ',num2str(R_polygenicscore)),'FitBoxToText','on');
xlabel('Average frequency positive allele'); ylabel('National IQ');
legend('show','Location','northeastoutside')
plot(X2,pfit,'color','k','DisplayName','regression');
hold off
if SNP_number > 28 || BARon == 'n'
    subplot(1,2,2);
else
    plotnum2=[((pl1+1)*4)+3,((pl1+1)*4)+4,((pl1+2)*4)+3,((pl1+2)*4)+4];
    subplot(pl,4,plotnum2);
end
hold on

for run23 = 2:6
    poprefnum=(popref.(cat(2,'pom',(num2str(run23)))))-1; color = col.(cat(2,'c',(num2str(run23))));
    plot(X3(poprefnum),Y3(poprefnum),'x','MarkerSize',12,'MarkerFaceColor',color,'MarkerEdgeColor',color,...
        'DisplayName',cell2mat(plotname(run23-1)));
end

ylim([60 110]); title({'';'Piffer ommisions (ITU,CDX,GIH)'});
annotation('textbox',[0.92 0 0.3 0.3],'String',cat(2,'r_2 (om) = ',num2str(R_polygenicscoreOM)),'FitBoxToText','on');
xlabel('Average frequency positive allele'); ylabel('National IQ');
legend('show','Location','northeastoutside')
plot(X3,pfit1,'color','k','DisplayName','regression');
disp '11 - Scatterplot and regression of polygenic score by sub population'; 
disp(cat(2,'r_1  for polygenic score vs IQ (all) = ',num2str(R_polygenicscore)));
disp(cat(2,'r_2 for polygenic score vs IQ (ommitted CDX, ITU, GIH) = ',num2str(R_polygenicscoreOM))); toc; disp ' ';

%% FA results w/ MLE estimates

figure('name','FA results w / MLE estimates')
set(gcf, 'Position', [50 10 1250 720])
X2=cell2mat(PolygenicScore(:,5)); Y2=cell2mat(PolygenicScore(:,4)); pcoef = polyfit(X2,Y2,1);
pfit = polyval(pcoef,X2); presid = Y2 - pfit; SSresid = sum(presid.^2); SStotal = (length(Y2)-1) * var(Y2);
rsq = 1 - SSresid/SStotal; R_fa = rsq^0.5; ax = subplot(1,2,1);
hold on
plot(ax,cell2mat(PolygenicScore(:,5)),cell2mat(PolygenicScore(:,4)),'x');
xlabel('FA Score'); ylabel('National IQ'); ylim([60 110]);
plot(cell2mat(PolygenicScore(:,5)),pfit,'color','k');
annotation('textbox',[0.91 0.40 .3 .3],'String',cat(2,'r_3 (all) = ',num2str(R_fa)),'FitBoxToText','on');
annotation('textbox',[0.47 0.55 .3 .3],'String',Fa_loading,'FitBoxToText','on');
annotation('textbox',[0.47 0.58 .3 .3],'String','FA Loading','FitBoxToText','on','EdgeColor','none');
X3=cell2mat(PolygenicScoreOM(:,5)); Y3=cell2mat(PolygenicScoreOM(:,4)); pcoef = polyfit(X3,Y3,1);
pfit1 = polyval(pcoef,X3); presid = Y3 - pfit1; SSresid = sum(presid.^2);
SStotal = (length(Y3)-1) * var(Y3); rsq = 1 - SSresid/SStotal; R_faOM = rsq^0.5; ax = subplot(1,2,2);
hold on
plot(ax,cell2mat(PolygenicScoreOM(:,5)),cell2mat(PolygenicScoreOM(:,4)),'x');
xlabel('FA Score'); ylabel('National IQ'); ylim([60 110]);
plot(cell2mat(PolygenicScoreOM(:,5)),pfit1,'color','k','DisplayName','regression');
annotation('textbox',[0.91 0.35 .3 .3],'String',cat(2,'r_4 (om) = ',num2str(R_faOM)),'FitBoxToText','on');
disp '12 - Scatterplot and regression of factor score by sub population';
disp(cat(2,'r_3  for factor score vs IQ (all) = ',num2str(R_fa)));
disp(cat(2,'r_4 for factor score vs IQ (ommitted CDX, ITU, GIH) = ',num2str(R_faOM))); toc; disp ' ';

%% Clear junk / placeholder variables

clear('A','aANC','aDER','run1','run2','run3','run4','run5','run6','run7','run8','run9',...
    'run10','run11','run12','run13','run14','run15','run16','run17','run18','run19','run20',...
    'run21','run22','run23','stat1','stat2','stat3','charlim','color','ax','col','dev','ef','impact',...
    'plotname','Anovacat','popref','prompt','statgroup','SSresid','SStotal','rsq','presid','Fa_psi',...
    'Fa_T','dim','err','X','X2','X3','Y','Y2','Y3','Ypos','YposFA','titleRS','superef','AA','REF',...
    'Tableout1','hBar','a21','a22','a23','a24','a25','Allele','hbar','impact','Known','pcoef','pfit','pfit1',...
    'pl','pl1','plotnum','plotnum2','poprefnum','ValueKnown','ValueKnown1','prevdir',...
    'Ab','Aref','Adir','Az','n','m','raw1','raw2','rsc','rsc1','BARon');

disp 'Completed'; toc;
