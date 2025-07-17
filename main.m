u1 = [2, 4, 0];
u2 = [5, 2, 0];
u3 = [7.98, 4, 0];
u4 = [6.5, 8, 0];
u5 = [5.5, 4, 0];
u6 = [-2, 1, 0];
u7 = [-3, 4, 0];
u8 = [-3, 1, 0];
Usuarios = {
      u1,
      u2,
      u3,
      u4,
      u5,
      u6,
      u7,
      u8
};
U = length(Usuarios);  % Número total de usuários
uhat = 6;              % Primeiro usuario no room2 (û)

VLC = [2.5, 2.5, 3]; # Meio do room1

espelhos = {};
refratores = {};

# Coloca as coordenadas de cada espelho e refrator ao redor do meio da parede entre rooms 1 e 2 (0, 2.5, 1.5)
for i = 1:5
  for j = 1:10
    y = 2.05 + (j - 1) * 0.1;
    z = 1.7 - (i - 1) * 0.1;

    if mod(i + j, 2) == 0
      espelhos{end+1} = [0, y, z];
    else
      refratores{end+1} = [0, y, z];
    endif
  endfor
endfor

# declarando do lado de fora, já que teremos que otimizar esses valores depois (Parte IV do artigo)
yaw = 45;
roll = 10;
Eta_c = 1.5;

lambda = 510; # muda entre 510 e 670 dependendo do teste
iota = 0; # 1 = com linha de visao; 0 = sem linha de visao
channel_gain = [];
Gamma = [];
for i = 1:U
  [channel_gain(i), Gamma(i)] = H(Usuarios{i}, VLC, espelhos, refratores, yaw, roll, Eta_c, lambda, iota);
endfor

#Ordena o primeiro grupo (1 até û-1)
[channel_gain1_ord, idx1] = sort(channel_gain(1:uhat-1), 'descend');
idx1_real = idx1;  #indices reais do primeiro grupo
#Ordena o segundo grupo (û até U)
[channel_gain2_ord, idx2] = sort(channel_gain(uhat:U), 'descend');
idx2_real = idx2 + (uhat - 1);  #Ajusta para o indice global
#Junta os dois grupos ordenados
channel_gain_ordenado = [channel_gain1_ord, channel_gain2_ord];
idx_total = [idx1_real, idx2_real];

channel_gain
Gamma
channel_gain_ordenado
idx_total
