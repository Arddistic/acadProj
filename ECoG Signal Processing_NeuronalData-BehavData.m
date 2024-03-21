clear all;
clc;

load('BME_526_HW2_NeuralData.mat');       % Loaded Neural data
load('BME_526_HW2_BehaviorData.mat');     % Loaded Behavior data

%% Data info

AcP_max = 0;
startTime = -1.0;
endTime = 2.5;
                                         % Working on Channel 20, 25 & 57

neuronCol = Channels.Chan20(:,2);        % Neuron values
uniqueNeuron = unique(neuronCol)         % finding unique neurons

Events = tm_times(1,:);                  % Setting an event

%% Event alignment of Action potentials

% Action potentials corresponding to neuron 1
% AcP_Val = Channels.Chan25(Channels.Chan25(:,2) == 1,3); 
% AcP_Val = Channels.Chan57(Channels.Chan57(:,2) == 1,3);
AcP_Val = Channels.Chan20(Channels.Chan20(:,2) == 1,3);


for k = 1:size(Events,2)
    
    Neuron1_aligned = AcP_Val - Events(k);
    
    AcP_Index = Neuron1_aligned >= startTime & Neuron1_aligned <= endTime;
    
    if sum(AcP_Index) > AcP_max
        AcP_max = sum(AcP_Index);
    end
    
    AcP_cell{k,:} = Neuron1_aligned(AcP_Index);
end

AcP_mtrx = zeros(size(Events,2), AcP_max);       % Cell array to matrix

for i=1:size(AcP_mtrx,1)
    
    AcP_mtrx(i,1:size(AcP_cell{i,:},1)) = AcP_cell{i,:};
end

AcP_mtrx(AcP_mtrx == 0) = "NaN";

%% Plots

figure;                                          % Raster plot
subplot(1, 2, 1)
plot(AcP_mtrx', 1:size(AcP_mtrx,1), 'b.')
xlim([startTime endTime])
title('Channel 20 N1 thumb-mid-finger event')
                                                 % PETH
subplot(1, 2, 2)
histogram(AcP_mtrx', 100, 'FaceColor', 'y','EdgeColor','y'  )
% changing histogram color (Colors make everything interesting!)
xlim([startTime endTime])
title('thumb-index-mid-finger')
line([startTime endTime],[th th], 'Color', 'm', 'LineStyle', '-.')

%% Electrode mapping 

load('ElectrodeMap');

channel = [20; 25; 57];
l = length(MapStruct.CerebusElectrode);

for d = 1:l;
    for c = 1:3;
        if MapStruct.CerebusChannel(d) == channel(c);
            E_cerebus(c) = MapStruct.CerebusElectrode(d);
        end

        if MapStruct.TDTChannel(d) == channel(c);
            E_TDT(c) = MapStruct.CerebusElectrode(d);
        end
    end
end

% For easy visualization of result

row = {'Channel', 'Electrode_cere', 'Electrode_TDT'};
y = [channel E_cerebus E_TDT];
x = array2table(y, 'VariableNames', row)