clc;
clear all;
close all;

% SYMBOLS
% R - RESISTANCE
% L - INDUCTANCE
% C - CAPACITANCE
% V - VOLTAGE SOURCE
% I - CURRENT SOURCE
% E - VCVS
% G - VCCS
% H - CCVS
% F - CCCS

prompt = 'Please enter the file name: ';
fileName = input(prompt,'s');
netlist_fileID = fopen(fileName);
netlist = textscan(netlist_fileID,'%s %s %s %s %s %s ');
fclose(netlist_fileID);
fileID1=fopen('Element_indep.txt','wt+'); %Create an empty text file Element_indep.txt for passive elements and independent sources
fileID2=fopen('VCVS.txt','wt+'); %Create an empty text file VCVS.txt for voltage controlled voltage sources
fileID3=fopen('VCCS.txt','wt+'); %Create an empty text file VCCS.txt for voltage controlled current sources
fileID4=fopen('CCVS.txt','wt+'); %Create an empty text file CCVS.txt for current controlled voltage sources
fileID5=fopen('CCCS.txt','wt+'); %Create an empty text file CCCS.txt for current controlled current sources

%---------------------------------------------------------------------------------------------------------------------------------------

msize = 0; % size of matrix G,X,C,dv/dt
num_R = 0; %number of resistances in circuit
num_C = 0; %number of capacitances in circuit
num_L = 0; %number of inductance in circuit 
num_V = 0; %number of independent voltage source in circuit
num_I = 0; %number of independent current source in circuit
num_Elements=0; %Number of passive elements
num_Nodes=0; %Number of nodes, excluding ground (node 0)
num_VCVS=0; %Number of voltage controlled voltage sources
num_VCCS=0; %Number of voltage controlled current sources
num_CCVS=0; %Number of current controlled voltage sources
num_CCCS=0; %Number of current controlled current sources

%------------------------------------------------------------------------------------------------------------------

for i=1:length(netlist{1})
    s=netlist{1}{i};
    display(s);
    switch(s(1))
        case{'R','L','C','V','I'} %For passive elements and independent sources
            fprintf(fileID1,[netlist{1}{i} ' ' netlist{2}{i} ' ' ...
                netlist{3}{i} ' ' netlist{4}{i} '\n']);
        case{'E'} %For voltage controlled voltage sources
            fprintf(fileID2,[netlist{1}{i} ' ' netlist{2}{i} ' ' ...
                netlist{3}{i} ' ' netlist{4}{i} ' ' netlist{5}{i} ' ' ...
                netlist{6}{i} '\n']);
        case{'G'} %For voltage controlled current sources
            fprintf(fileID3,[netlist{1}{i} ' ' netlist{2}{i} ' ' ...
                netlist{3}{i} ' ' netlist{4}{i} ' ' netlist{5}{i} ' ' ...
                netlist{6}{i} '\n']);
        case{'H'} %For current controlled voltage sources
            fprintf(fileID4,[netlist{1}{i} ' ' netlist{2}{i} ' ' ...
                netlist{3}{i} ' ' netlist{4}{i} ' ' netlist{5}{i} '\n']);
        case{'F'} %For current controlled current sources
            fprintf(fileID5,[netlist{1}{i} ' ' netlist{2}{i} ' ' ...
                netlist{3}{i} ' ' netlist{4}{i} ' ' netlist{5}{i} '\n']);
    end
end

%----------------------------------------------------------------------------------------------------------------------

[Name,N1,N2,value]=textread('Element_indep.txt','%s %s %s %s');
for i=1:length(Name)
    switch(Name{i}(1))
        case{'R','L','C'}
            num_Elements=num_Elements+1;
            Element(num_Elements).Name=Name{i};
            Element(num_Elements).Node1=str2num(N1{i});
            Element(num_Elements).Node2=str2num(N2{i});
            Element(num_Elements).Value=str2double(value{i});
            if(Name{i}(1)=='L')
                num_L=num_L+1;
                Inductor(num_L).Name=Name{i};
                Inductor(num_L).N1=str2num(N1{i});
                Inductor(num_L).N2=str2num(N2{i});
                Inductor(num_L).Value=str2double(value{i});
            elseif(Name{i}(1)=='C')
                num_C=num_C+1;
                Capacitor(num_C).Name=Name{i};
                Capacitor(num_C).N1=str2num(N1{i});
                Capacitor(num_C).N2=str2num(N2{i});
                Capacitor(num_C).Value=str2double(value{i});
            else
                num_R=num_R+1;
                Resistance(num_R).Name=Name{i};
                Resistance(num_R).N1=str2num(N1{i});
                Resistance(num_R).N2=str2num(N2{i});
                Resistance(num_R).Value=str2double(value{i});
            end
        case{'V'}
            num_V=num_V+1;
            Volt_source(num_V).Name=Name{i};
            Volt_source(num_V).Node1=str2num(N1{i});
            Volt_source(num_V).Node2=str2num(N2{i});
            Volt_source(num_V).Value=str2double(value{i});
            
        case{'I'}
            num_I=num_I+1;
            Current_source(num_I).Name=Name{i};
            Current_source(num_I).Node1=str2num(N1{i});
            Current_source(num_I).Node2=str2num(N2{i});
            Current_source(num_I).Value=str2double(value{i});
    end
    num_Nodes=max(str2num(N1{i}),max(str2num(N2{i}),num_Nodes));
end

%--------------------------------------------------------------------------------------------------------------
%VCVS
[Name1,N1,N2,NC1,NC2,Gain]=textread('VCVS.txt','%s %s %s %s %s %s');
num_VCVS=length(Name1);
for i=1:num_VCVS
    VCVS(i).Name=Name1{i};
    VCVS(i).N1=str2num(N1{i});
    VCVS(i).N2=str2num(N2{i});
    VCVS(i).NC1=str2num(NC1{i});
    VCVS(i).NC2=str2num(NC2{i});
    VCVS(i).Gain=str2double(Gain{i});
    num_Nodes=max(str2num(N1{i}),max(str2num(N2{i}),num_Nodes));
end
%--------------------------------------------------------------------------------------------------------------
%VCCS
[Name2,N1,N2,NC1,NC2,Transconductance]=textread('VCCS.txt','%s %s %s %s %s %s');
num_VCCS=length(Name2);
for i=1:num_VCCS
    VCCS(i).Name=Name2{i};
    VCCS(i).N1=str2num(N1{i});
    VCCS(i).N2=str2num(N2{i});
    VCCS(i).NC1=str2num(NC1{i});
    VCCS(i).NC2=str2num(NC2{i});
    VCCS(i).Transconductance=str2double(Transconductance{i});
    num_Nodes=max(str2num(N1{i}),max(str2num(N2{i}),num_Nodes));
end
%--------------------------------------------------------------------------------------------------------------
%CCVS
[Name3,N1,N2,Vcontrol,Transresistance]=textread('CCVS.txt','%s %s %s %s %s');
num_CCVS=length(Name3);
for i=1:num_CCVS
    CCVS(i).Name=Name3{i};
    CCVS(i).N1=str2num(N1{i});
    CCVS(i).N2=str2num(N2{i});
    CCVS(i).Vcontrol=Vcontrol{i};
    CCVS(i).Transresistance=str2double(Transresistance{i});
    num_Nodes=max(str2num(N1{i}),max(str2num(N2{i}),num_Nodes));
end
%--------------------------------------------------------------------------------------------------------------------
%CCCS
[Name4,N1,N2,Vcontrol,Gain]=textread('CCCS.txt','%s %s %s %s %s');
num_CCCS=length(Name4);
for i=1:num_CCCS
    CCCS(i).Name=Name4{i};
    CCCS(i).N1=str2num(N1{i});
    CCCS(i).N2=str2num(N2{i});
    CCCS(i).Vcontrol=Vcontrol{i};
    CCCS(i).Gain=str2double(Gain{i});
    num_Nodes=max(str2num(N1{i}),max(str2num(N2{i}),num_Nodes));
end
%---------------------------------------------------------------------------------------------------------------------------
%Matrix G, C size.
msize = num_Nodes + num_L + num_VCVS + num_V + num_VCVS + num_CCVS + num_CCVS + num_CCCS;
%----------------------------------------------------------------------------------------------------------------------------
%Initializing G, X, and C Matrix
G = zeros(msize);
C = zeros(msize);
%X = zeros(msize,2);
B = zeros(msize,1);
goa = zeros(msize,1);

a = num_Nodes;
mah = 1;
%----------------------------------------------------------------------------------------------------------------------------
%Populating G and C Matrix 
for i=1:num_Elements
    switch(Element(i).Name(1))
        case 'R'
            if(Element(i).Node1 == 0 || Element(i).Node2 == 0)
                if(Element(i).Node1 == 0)
                    G(Element(i).Node2,Element(i).Node2) = G(Element(i).Node2,Element(i).Node2) + (1/(Element(i).Value));
                    X{Element(i).Node2,1} =strcat('V_',num2str(Element(i).Node2));
                else
                    G(Element(i).Node1,Element(i).Node1) = G(Element(i).Node1,Element(i).Node1) + (1/(Element(i).Value));
                    X{Element(i).Node1,1} =strcat('V_',num2str(Element(i).Node1));
                end
            else
                G(Element(i).Node1,Element(i).Node1) = G(Element(i).Node1,Element(i).Node1) + (1/(Element(i).Value));
                G(Element(i).Node1,Element(i).Node2) = G(Element(i).Node1,Element(i).Node2) - (1/(Element(i).Value));
                X{Element(i).Node1,1} =strcat('V_',num2str(Element(i).Node1));
                G(Element(i).Node2,Element(i).Node1) = G(Element(i).Node2,Element(i).Node1) - (1/(Element(i).Value));
                G(Element(i).Node2,Element(i).Node2) = G(Element(i).Node2,Element(i).Node2) + (1/(Element(i).Value));
                X{Element(i).Node2,1} =strcat('V_',num2str(Element(i).Node2));
            end
        case 'C'
            if(Element(i).Node1 == 0 || Element(i).Node2 == 0)
                if(Element(i).Node1 == 0)
                    C(Element(i).Node2,Element(i).Node2) = C(Element(i).Node2,Element(i).Node2) + (Element(i).Value);
                    X{Element(i).Node2,1} =strcat('V_',num2str(Element(i).Node2));
                else
                    C(Element(i).Node1,Element(i).Node1) = C(Element(i).Node1,Element(i).Node1) + (Element(i).Value);
                    X{Element(i).Node1,1} =strcat('V_',num2str(Element(i).Node1));
                end
            else
                C(Element(i).Node1,Element(i).Node1) = C(Element(i).Node1,Element(i).Node1) + (Element(i).Value);
                C(Element(i).Node1,Element(i).Node2) = C(Element(i).Node1,Element(i).Node2) - (Element(i).Value);
                C(Element(i).Node2,Element(i).Node1) = C(Element(i).Node2,Element(i).Node1) - (Element(i).Value);
                C(Element(i).Node2,Element(i).Node2) = C(Element(i).Node2,Element(i).Node2) + (Element(i).Value);
                X{Element(i).Node2,1} =strcat('V_',num2str(Element(i).Node2));
                X{Element(i).Node1,1} =strcat('V_',num2str(Element(i).Node1));
            end
         case 'L'
                 if(Element(i).Node1 == 0 || Element(i).Node2 == 0)
                     a = a + 1;
                     if(Element(i).Node2 == 0)
                         G(Element(i).Node1,a) = 1;
                         G(a,Element(i).Node1) = 1;
                         C(a,a) = -(Element(i).Value);
                         X{a,1} =strcat('I_L',num2str(mah));
                         mah = mah + 1;
                     else
                         G(Element(i).Node2,a) = -1;
                         G(a,Element(i).Node2) = -1;
                         C(a,a) = -(Element(i).Value);
                         X{a,1} =strcat('I_L',num2str(mah));
                         mah = mah + 1;
                     end
                 else
                     a = a + 1;
                     G(Element(i).Node1,a) = 1;
                     G(Element(i).Node2,a) = -1;
                     G(a,Element(i).Node1) = 1;
                     G(a,Element(i).Node2) = -1;
                     C(a,a) = -(Element(i).Value);
                     X{a,1} =strcat('I_L',num2str(mah));
                     mah = mah + 1;
                 end
    end
end
%--------------------------------------------------------------------------------------------------------------------------------------
%Voltage Source populating
mah = 1;
for i = 1:num_V
    a = a + 1;
    if(Volt_source(i).Node1 == 0 || Volt_source(i).Node2 == 0)
        if(Volt_source(i).Node1 == 0)
            G(Volt_source(i).Node2,a) = 1;
            G(a,Volt_source(i).Node2) = -1;
            B(a) = Volt_source(i).Value;
            X{a,1} =strcat('I_Vs',num2str(mah));
            mah = mah + 1;
        else
            G(Volt_source(i).Node1,a) = -1;
            G(a,Volt_source(i).Node1) = 1;
            B(a) = Volt_source(i).Value;
            X{a,1} =strcat('I_Vs',num2str(mah));
            mah = mah + 1;
        end
    else
        B(a) = Volt_source(i).Value;
        G(Volt_source(i).Node1,a) = -1;
        G(Volt_source(i).Node2,a) = 1;
        G(a,Volt_source(i).Node1) = 1;
        G(a,Volt_source(i).Node2) = -1;
        X{a,1} =strcat('I_Vs',num2str(mah));
        mah = mah + 1;
    end
end
%--------------------------------------------------------------------------------------------------------------------------------------
%Current Source populating
for i = 1:num_I
    if(Current_source(i).Node1 == 0 || Current_source(i).Node2 == 0)
        if(Current_source(i).Node1 == 0)
            B(Current_source(i).Node2) = B(Current_source(i).Node2) - (Current_source(i).Value);
        else
            B(Current_source(i).Node1) = B(Current_source(i).Node1) + Current_source(i).Value;
        end
    else
        B(Current_source(i).Node1) = B(Current_source(i).Node1) - (Current_source(i).Value);
        B(Current_source(i).Node2) = B(Current_source(i).Node2) + Current_source(i).Value;
    end
end
%--------------------------------------------------------------------------------------------------------------------------------------
%VCVS populating
mah = 1;
for i=1:num_VCVS
    a=a+1;
    if(VCVS(i).N1 == 0 || VCVS(i).N2 == 0)
        if(VCVS(i).N2 == 0)
            G(VCVS(i).N1,a)= -1;
            G(a,VCVS(i).N1)= 1;
            G(a,VCVS(i).NC1)= -(VCVS(i).Gain);
            G(a,VCVS(i).NC2)= (VCVS(i).Gain);
            B(a)=0;
            X{a,1} =strcat('I_VCVS',num2str(mah));
            mah = mah + 1;
        else
            G(VCVS(i).N2,a)= 1;
            G(a,VCVS(i).N2)= -1;
            G(a,VCVS(i).NC1)= -(VCVS(i).Gain);
            G(a,VCVS(i).NC2)= (VCVS(i).Gain);
            B(a)=0;
            X{a,1} =strcat('I_VCVS',num2str(mah));
            mah = mah + 1;
        end
    elseif((VCVS(i).N2 == 0 && VCVS(i).NC2 == 0)  || (VCVS(i).N1 == 0 && VCVS(i).NC1 == 0))
        if(VCVS(i).N2 == 0 && VCVS(i).NC2 == 0)
            G(VCVS(i).N1,a)= -1;
            G(a,VCVS(i).N1)= 1;
            G(a,VCVS(i).NC1)= -(VCVS(i).Gain);
            B(a)=0;
            X{a,1} =strcat('I_VCVS',num2str(mah));
            mah = mah + 1;
        else
            G(VCVS(i).N2,a)= 1;
            G(a,VCVS(i).N2)= -1;
            G(a,VCVS(i).NC2)= (VCVS(i).Gain);
            B(a)=0;
        end
    else
         G(VCVS(i).N1,a)= -1;
         G(VCVS(i).N2,a)= 1;
         G(a,VCVS(i).NC1)= -(VCVS(i).Gain);
         G(a,VCVS(i).N1)= 1;
         G(a,VCVS(i).N2)= -1;
         G(a,VCVS(i).NC2)= (VCVS(i).Gain);
         B(a)=0;
         X{a,1} =strcat('I_VCVS',num2str(mah));
         mah = mah + 1;
    end
    
end
    
%-------------------------------------------------------------------------------------------------------------------------------------
%VCCS populating

for i=1:num_VCCS
    if(VCCS(i).N1 == 0 || VCCS(i).N2 == 0)
        if(VCCS(i).N2 == 0)
            G(VCCS(i).N1,VCCS(i).NC1)= -(VCCS(i).Transconductance);
            G(VCCS(i).N1,VCCS(i).NC2)= (VCCS(i).Transconductance);
            X{Element(i).Node1,1} =strcat('V_VCCS',num2str(Element(i).Node1));     
        else
            G(VCCS(i).N2,VCCS(i).NC1)= (VCCS(i).Transconductance);
            G(VCCS(i).N2,VCCS(i).NC2)= -(VCCS(i).Transconductance);
            X{Element(i).Node2,1} =strcat('V_VCCS',num2str(Element(i).Node2));
        end
    elseif((VCCS(i).N2 == 0 && VCCS(i).NC2 == 0)  || (VCCS(i).N1 == 0 && VCCS(i).NC1 == 0))
        if(VCCS(i).N2 == 0 && VCCS(i).NC2 == 0)
            G(VCCS(i).N1,VCCS(i).NC1)= -(VCCS(i).Transconductance);
            X{Element(i).Node1,1} =strcat('V_VCCS',num2str(Element(i).Node1)); 
        else
            G(VCCS(i).N2,VCCS(i).NC2)= -(VCCS(i).Transconductance);
            X{Element(i).Node2,1} =strcat('V_VCCS',num2str(Element(i).Node2));
        end
    else
        G(VCCS(i).N1,VCCS(i).NC1)= -(VCCS(i).Transconductance);
        G(VCCS(i).N2,VCCS(i).NC1)= VCCS(i).Transconductance;
        G(VCCS(i).N1,VCCS(i).NC2)= VCCS(i).Transconductance;
        G(VCCS(i).N2,VCCS(i).NC2)= -(VCCS(i).Transconductance);
        X{Element(i).Node1,1} =strcat('V_VCCS',num2str(Element(i).Node1));
        X{Element(i).Node2,1} =strcat('V_VCCS',num2str(Element(i).Node2));
    end
    
end

%--------------------------------------------------------------------------------------------------------------------------------------
%CCVS populating
mah = 1;
for i=1:num_CCVS
    a=a+1;
    b=a;
    a=b+1; %updating value of a
    G(CCVS(i).N1,b)=1;
    G(CCVS(i).N2,b)=-1;
    G(CCVS(i).NC1,a)=1;
    G(CCVS(i).NC2,a)=-1;
    G(b,CCVS(i).N1)=1;
    G(b,CCVS(i).N2)=-1;
    G(a,CCVS(i).NC1)=1;
    G(a,CCVS(i).NC2)=-1;
    G(b,a)=-CCVS(i).Transresistance;
    X{a,1} =strcat('Iout_CCVS',num2str(mah));
    X{b,1} =strcat('Iin_CCVS',num2str(mah));
    mah = mah + 1;
end
%--------------------------------------------------------------------------------------------------------------------------------------
%CCCS populating
mah = 1;
for i=1:num_CCCS
    a=a+1;
    G(CCCS(i).N1,a)= CCCS(i).Gain;
    G(CCCS(i).N2,a)= -(CCCS(i).Gain);
    G(CCCS(i).NC1,a)= 1;
    G(CCCS(i).NC2,a)= -1;
    G(a,CCCS(i).NC1)= 1;
    G(a,CCCS(i).NC2)= -1;
    X{a,1} =strcat('I_CCCS',num2str(mah));
    mah = mah + 1;
end
%---------------------------------------------------------------------------------------------------------------------------------
%Solving the Matrices

% B = G*X + C*dX;
% C*dX = B - G*X;  ------ Eq 1
% Applying Backward Euler Formula
% C*X_(n+1)= C*X_n + h*C*dX(n+1) ----------- Eq 2
% Substituting Eq 1 into Eq 2
% (C + h*G)*X_(n+1) = C*X_n + h*B_(n+1)

h = input('Enter Value of time step (h):');
pun = C+h*G;
mum = inv(pun);
goa =  goa + (mum*(C*goa + h*B));


