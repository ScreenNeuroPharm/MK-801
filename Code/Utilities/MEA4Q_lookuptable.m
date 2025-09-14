function [MEA4Q] = MEA4Q_lookuptable()
%UNTITLED Summary of this function goes here
%   Detailed explanation goes here

    clear all
    clc
    MEA4Q = string(ones(60,2));
    letter = 'A':'O';
    letter(letter=='J' | letter =='N')= [];
    count = 1;
    
    for i = 1:length(letter)
        if letter(i) == 'A' | letter(i) == 'O'
            for k = [2,3,5,6]
                MEA4Q(count,1) = strcat(letter(i),'0',char(string(k)));
                count = count+1;
            end
        elseif letter(i) == 'B' | letter(i) == 'M'
            for k = [1:3,5:7]
                MEA4Q(count,1) = strcat(letter(i),'0',char(string(k)));
                count = count+1;
            end
        elseif letter(i) == 'E' | letter(i) == 'I'
             for k = [1,2,6,7]
                MEA4Q(count,1) = strcat(letter(i),'0',char(string(k)));
                count = count+1;
             end
        elseif letter(i) == 'F' | letter(i) == 'G' 
            k = 4;
            MEA4Q(count,1) = strcat(letter(i),'0',char(string(k)));
            count = count+1;
        elseif letter(i) == 'H'
            for k = [1,4]
                MEA4Q(count,1) = strcat(letter(i),'0',char(string(k)));
                count = count+1;
             end
        else 
            for k = [1:7]
                MEA4Q(count,1) = strcat(letter(i),'0',char(string(k)));
                count = count+1;
            end
        end
    end
    
    mcs = [32 31 61 62 33 21 44 54 71 63 12 22 43 52 53 72 82 13 23 41 42 51 ...
            73 83 24 34 64 74 84 14 15 85 25 35 65 75 16 26 48 57 58 76 86 17 27 46 47 56 ...
            77 87 36 28 45 55 78 66 37 38 68 67]';
    MEA4Q(:,2)=mcs;
    
            


