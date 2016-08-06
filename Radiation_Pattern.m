function Radiation_Pattern()

%*****************************************
%************* CREATE FOLDER *************
%*****************************************
%create folder to save Diagrams and Infos
currentDate = date;
myFolder = 'Radiation Pattern _ ';
newFolder = [myFolder currentDate];
mkdir(newFolder)
cd(newFolder)
%write read me.txt in newFolder
readMe = fopen('READ ME.txt', 'wt');
fprintf(readMe,'This folder, is created to store radiation diagrams!\n');
fprintf(readMe,'Please not write here\n');
fprintf(readMe,'Thank you!');
fclose(readMe);


