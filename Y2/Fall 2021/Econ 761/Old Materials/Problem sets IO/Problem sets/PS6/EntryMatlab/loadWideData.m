
clear all;

numAirlines=31;
numCarVar=9;
load airlines5test;
dlength=length(data(:,1));

pax=zeros(dlength,numAirlines);
pot96=zeros(dlength,numAirlines);
pot97=zeros(dlength,numAirlines);
enter=zeros(dlength,numAirlines);
city2=zeros(dlength,numAirlines);
enter96=zeros(dlength,numAirlines);
numroute=zeros(dlength,numAirlines);
sharepax=zeros(dlength,numAirlines);
sharepaxdist=zeros(dlength,numAirlines);   
   

i=2;
j=1;
while i<(numAirlines*numCarVar+1),
   pax(:,j)=data(:,i);
   i=i+1;
   pot96(:,j)=data(:,i);
   i=i+1;
	pot97(:,j)=data(:,i);
   i=i+1;
   enter(:,j)=data(:,i);
   i=i+1;
   city2(:,j)=data(:,i);
   i=i+1;
   enter96(:,j)=data(:,i);
   i=i+1;
   numroute(:,j)=data(:,i);
   i=i+1;
   sharepax(:,j)=data(:,i);
   i=i+1;
   sharepaxdist(:,j)=data(:,i);
   i=i+1;
   j=j+1;
end;


citypair=data(:,1);
city972=data(:,numCarVar*numAirlines+2);
distance=data(:,numCarVar*numAirlines+3);
tourist=data(:,numCarVar*numAirlines+4);
basepop=data(:,numCarVar*numAirlines+5);
refpop=data(:,numCarVar*numAirlines+6);
paxtot=data(:,numCarVar*numAirlines+7);
pop=data(:,numCarVar*numAirlines+8);
totenter=data(:,numCarVar*numAirlines+9);
cityN2=data(:,numCarVar*numAirlines+10);
cityN1=data(:,numCarVar*numAirlines+11);
totpotential=data(:,numCarVar*numAirlines+12);
herfCityPair=data(:,numCarVar*numAirlines+13);
incumbents=data(:,numCarVar*numAirlines+14);
totsinglepot=data(:,numCarVar*numAirlines+15);
dist2=data(:,numCarVar*numAirlines+16);




%citypair "pax_1" "pot96_1" "pot97_1" "enter_1" "city2_1" "enter96_1" "numroute_1" "sharepax_1" "sharepaxdist_1" "city972" "base_cty" "ref_cty" "distance" "tourist" "base_apt" "ref_apt" "basepop" "refpop" "paxtot" "citypairpop" "totenter" "cityN2" "cityN1" "dist2"





