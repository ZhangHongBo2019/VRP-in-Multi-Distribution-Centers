clc;
clear;
%% ���������ĵĳ�����������
%��������
load data.mat
%����λ�þ���
m=size(X,1);
D=zeros(m,m);
for i=1:m
    for j=1:m
        D(i,j)=norm(X(i,:)-X(j,:));
        D(j,i)=D(i,j);
        D(i,i)=eps;
    end
end
%�����������ĵ�λ��
nume=zeros(20,20);
for i=1:20
    for j=1:20 
        temp=[D(:,i) D(:,j)];
        [num,pp]=min(temp,[],2);
        nume(i,j)=sum(num);
    end
end
[mum1,index1]=min(nume);
[mun2,index2]=min(mum1,[],2);
w1=index1(index2);
w2=index2;
%�������������ķ���˿ͼ�
H=[];
S=[];
for i=1:20
    if D(i,w1)<D(i,w2)
        H=[H i];
    else 
        S=[S i];
    end
end
n2=size(S,2);
%��Ⱥ�㷨����С�ĳ������г�
%���ò���
Pop=30;%��Ⱥ��Ŀ
Alpha=1;%��Ҫ��ϵ��
Beta=1;%Beta���ܼ���ϵ��
gama=2;
Rho=0.15;%�ӷ���ϵ��
MAXGEN=30;%��������
Q=15;%��Ϣ���²���
W=9;%W:����������
V=10;%�������ݻ�
w=[2.5,1.5,1.8,2.0,0.8,1.5,1.0,2.5,3.0,1.7,0.6,0.2,2.4,1.9,2.0,0.7,0.5,2.2,3.1,0.1];%ÿ���ͻ�����Ļ�������
t=[1.5,3.8,0.5,3,2.6,3.6,1.4,2.4,2,3.4,2,1.2,0.5,0.8,1.3,1.6,1.7,0.5,0.8,1.4];%ÿ���ͻ�����Ļ�����ݻ�
load_w_S=0;
load_t_S=0;
Eta=1./D;%�������ӣ���Ϊ����ĵ���
Tau_S=ones(m,m);
Tabu_S=zeros(Pop,n2+10);
iter=1;
G_best_route_S=[MAXGEN,n2+10];
G_best_length_S=zeros(MAXGEN,1);
length_ave_S=zeros(MAXGEN,1);%����·�ߵ�ƽ������
%��ʼ���е���
while iter<=MAXGEN
    Tabu_S(:,1)=w2*ones(Pop,1);
 for i=1:Pop
        visited_S=Tabu_S(i,:);
        visited_S=visited_S(visited_S>0);
        to_visit_S=setdiff(S,visited_S);
        d=1;
 %��w2Ϊ��������
 while d<=n2
 if ~isempty(to_visit_S)
%���չ���ѡ��һ�����������ǻص��ֿ�
  x=to_visit_S;
 for k=1:length(to_visit_S)
  x(k)=(Tau_S(visited_S(end),to_visit_S(k))^Alpha)*(Eta(visited_S(end),to_visit_S(k))^Beta);
 end
 x=x/(sum(x)); 
%������ԭ��ѡȡ��һ������ 
XC=cumsum(x); 
Select=find(XC>=rand);
if isempty(Select)
    Select=w2;
    load_w_S=load_w_S+w(Select);
    load_t_S=load_w_S+t(Select);
else
load_w_S=load_w_S+w(to_visit_S(Select(1)));
load_t_S=load_t_S+t(to_visit_S(Select(1)));
end
c1=min((load_w_S)-9,0);c2=min((load_t_S)-10,0);
if c1<0&&c2<0
    Tabu_S(i,length(visited_S)+1)=to_visit_S(Select(1)); 
else
    Select=w1;
       j=j-1;
    load_w_S=0;
    load_t_S=0;
    Tabu_S(i,length(visited_S)+1)=Select;
end
 end
    visited_S=Tabu_S(i,:);
    visited_S=visited_S(visited_S>0);
    to_visit_S=setdiff(S,visited_S);
     if visited_S(end)~=w2
   Tabu_S(i,1:(length(visited_S)+1))=[visited_S,w2];
     end
           d=d+1;
 end
    load_w_S=0;
    x=[];   
 end
%% ���Ĳ���¼�������ֲ����������ֻ���ϵ�·��
L_S=zeros(Pop,1);
for i=1:Pop  
MM_S=Tabu_S(i,:);
R_S=MM_S(MM_S>0);
for j=1:(length(R_S)-1)
L_S(i)=L_S(i)+D(R_S(j),R_S(j+1));
end 
end
%������̾�������·�� 
[min_Length_S,index_S]=min(L_S);
G_best_length_S(iter)=min_Length_S;
G_best_route_S(iter,1:length(Tabu_S(index_S(1),:)))=Tabu_S(index_S(1),:);
length_ave_S(iter)=mean(L_S); 
disp(['��',num2str(iter),'��']);
iter=iter+1;
%% ���岽������Ϣ��
Delta_Tau_S=zeros(m,m); 
for i=1:Pop 
MM=Tabu_S(i,:); 
R=MM(MM>0);
for j=1:(length(R)-1) 
Delta_Tau_S(R(j),R(j+1))=Delta_Tau_S(R(j),R(j+1))+Q/L_S(i); 
end 
end 
Tau_S=(1-Rho).*Tau_S+Delta_Tau_S;
%% �����������ɱ����� 
Tabu_S=zeros(Pop,n2); 
load_w_S=0;
end
%% ���߲��������� 
[best_length_S,index]=min(G_best_length_S); 
best_route_S=G_best_route_S(index(1),:);
best_route_S=best_route_S(best_route_S>0);
disp('~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~');
disp(['���·��Ϊ��',num2str(best_route_S)]);
disp(['���·��Ϊ��',num2str(best_length_S)]);
%% �ڰ˲�������ɢ��ͼ��Ѳ�ι���ͼ
%����ɢ��ͼ������ע�������ĵ�λ��
figure(1)
plot(X(:,1),X(:,2),'o');
hold on
plot(X(best_route_S,1),X(best_route_S,2),'o-');
hold on
plot([X(w1,1),X(w2,1)],[X(w1,2),X(w2,2)],'rp','MarkerSize',9);
hold on
for i=1:n2
    text(X(S(i),1),X(S(i),2),['  ' num2str(S(i))]);
end
figure(2)
plot(1:MAXGEN,G_best_length_S) ;
hold on 
plot(1:MAXGEN,length_ave_S);
legend(' G_best_length_S ',' length_ave_S ');
