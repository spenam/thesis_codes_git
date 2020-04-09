using OrdinaryDiffEq, ForwardDiff, Distributions, RandomNumbers, NPZ

function initial_sampling(Eb,w1,mb)
    nums=rand(1)*2*pi
    x1=cos(nums[1])*sqrt(Eb/(0.5*mb*w1*w1))
    y1=sin(nums[1])*sqrt(2*mb*Eb)
    return x1,y1
end

function bolz_energy(a)
    return rand(Truncated(Exponential(1), 0., 10))
end
function bolz_freq(a)
    return rand()*0.4324555320336759+0.63
end
N=0

initial_position = [0.1, 0.,0.,0.,0.,0.,0.,0.,0.]
initial_momentum=[0. ,0.1,0.,0.,0.,0.,0.,0.,0.]
numeros=[2,4,6,8,10,16,20,26,30,36,40,50,60,80,100]
#numeros=[4,6]
promedios=100
solu=0
w=0
Es=0
count=zeros(size(numeros)[1])
for ii in 1:size(numeros)[1]
    println(ii)
    ii=Int(ii)
        
for jj in 1:promedios
    jj=Int(jj)
N=numeros[ii]
initial_position=ones(N+1)
initial_momentum=ones(N+1)

a=0.2
b=0.01

m=[1.,0.1,0.1,0.1,0.1,0.1,0.1,0.1,0.1]
w=[0.,0.8,0.8,0.8,0.8,0.8,0.8,0.8,0.8]
g=[0.,0.1,0.1,0.1,0.1,0.1,0.1,0.1,0.1]

X0=17.    
P0=1.

golden=(1+sqrt(big(5)))/2
Es=zeros(N+1)

time=10000.
#time=1.
tspan=(0., time)
#fig=figure()
w=zeros(N+1)
g=zeros(N+1)
E=zeros(N+1)
m=ones(N+1) .*0.1
m[1]=1.
w=w .* 0. .+ 1.
w=(bolz_freq.(w))
Es=(bolz_energy.(Es))



#println(w)
W=maximum(w)
g=g .* 0. .+ 0.3/N

dts=0.01*2. *pi/W
#Hsys(q,p) = p[1]^2. /(2. *m[1]) +m[1]*w[1]^2. *q[1]^2. /2. +0.5* -q[1]*g[2]*q[2]+p[2]^2. /(2. *m[2])+m[2]*w[2]^2. *q[2]^2. /2. + q[1]^2. *g[2]^2. /(2. *m[2]*w[2]^2.) -q[1]*g[3]*q[3]+p[3]^2. /(2. *m[3])+m[3]*w[3]^2. *q[3]^2. /2. + q[1]^2. *g[3]^2. /(2. *m[3]*w[3]^2.) -q[1]*g[4]*q[4]+p[4]^2. /(2. *m[4])+m[4]*w[4]^2. *q[4]^2. /2. + q[1]^2. *g[4]^2. /(2. *m[4]*w[4]^2.) -q[1]*g[5]*q[5]+p[5]^2. /(2. *m[5])+m[5]*w[5]^2. *q[5]^2. /2. + q[1]^2. *g[5]^2. /(2. *m[5]*w[5]^2.) -q[1]*g[6]*q[6]+p[6]^2. /(2. *m[6])+m[6]*w[6]^2. *q[6]^2. /2. + q[1]^2. *g[6]^2. /(2. *m[6]*w[6]^2.) -q[1]*g[7]*q[7]+p[7]^2. /(2. *m[7])+m[7]*w[7]^2. *q[7]^2. /2. + q[1]^2. *g[7]^2. /(2. *m[7]*w[7]^2.) -q[1]*g[8]*q[8]+p[8]^2. /(2. *m[8])+m[8]*w[8]^2. *q[8]^2. /2. + q[1]^2. *g[8]^2. /(2. *m[8]*w[8]^2.)-q[1]*g[9]*q[9]+p[9]^2. /(2. *m[9])+m[9]*w[9]^2. *q[9]^2. /2. + q[1]^2. *g[9]^2. /(2. *m[9]*w[9]^2.) 
function Hsys(q,p,N0)
    result=p[1]^2. /(2. *m[1]) - a*q[1]^2. /2. +b*q[1]^4. /4.
    for i in 1:N0
        result +=  q[1]*g[i+1]*q[i+1]+p[i+1]^2. /(2. *m[2])+m[i+1]*w[i+1]^2. *q[i+1]^2. /2.
    end
    return result
end
pdot(dp,p,q,params,t) = ForwardDiff.gradient!(dp, q->-Hsys(q, p,N), q)
qdot(dq,p,q,params,t) = ForwardDiff.gradient!(dq, p-> Hsys(q, p,N), p)



initial_position=initial_position .* 0.
initial_momentum=initial_momentum .* 0.
initial_position[1]=0.0 #xz0[jj]
initial_momentum[1]=0.0 #pz1[jj]

for i in 1:N
    i=Int(i)
    initial_position[i], initial_momentum[i]=initial_sampling(Es[i],w[i],m[i])
end

 
prob = DynamicalODEProblem(pdot, qdot, initial_momentum, initial_position, tspan)
@time sol = solve(prob, CalvoSanz4(), dt=dts)
    solu=sol

for j in 1:size(sol[N+2,:])[1]-1
    if sol[N+2,j]<0 && sol[N+2,j+1]>0
        count[ii]=count[ii]+1
    elseif sol[N+2,j]>0 && sol[N+2,j+1]<0
        count[ii]=count[ii]+1
    end
end
end
count=count ./promedios
end
println(count)
npzwrite("numeros_sin_extra.npy", numeros)
npzwrite("counts_sin_extra.npy",counts)
