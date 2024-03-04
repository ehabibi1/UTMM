%MATLAB R2021a

addpath('C:\Users\ehabibi\Downloads\fiji-win64\Fiji.app\scripts');
ImageJ

IJM.getDatasetAs('ecad');
IJM.show('ecad')
ecad_8b = uint8(ecad); % convert to uint8
IJM.show('ecad_8b')
