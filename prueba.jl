#julia --threads 4
using Base.Threads
using SparseArrays
using LoopVectorization
#print(Threads.nthreads())
#print("\n")
#print(Threads.threadid())
#print("\n")
K=spzeros(Threads.nthreads()+1,Threads.nthreads()+1)
print("This program is using ", Threads.nthreads(), " threads","\n")
#my_lock = ReentrantLock();
#lock(my_lock) 
@sync Threads.@threads for i = 1:Threads.nthreads()
   # print("i = $i on thread $(Threads.threadid()) \n")
    Aelem=[Threads.threadid() 0; 0 Threads.threadid()]
    K[Threads.threadid():Threads.threadid()+1,Threads.threadid():Threads.threadid()+1]=Aelem
end
#unlock(my_lock)
print(K)
