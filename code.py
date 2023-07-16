from operator import truediv
from re import T
import numpy as np
import random
from queue import PriorityQueue
import networkx as nx
from queue import Queue

#Gaussian Graph Generator
def Gaussian(m):

    G = nx.DiGraph()

    num = 1

    for k in range(1,m):

        for i in range(k+1,m+1):

            node1 = "T{},{}".format(k,k)

            node2 = "T{},{}".format(k,i)

            G.add_edge(node1,node2)

        if not k+1 == m:

            for j in range(k+1, m+1):

                t = k+1

                node1 = "T{},{}".format(k,j)

                node2 = "T{},{}".format(t,j)

                G.add_edge(node1, node2)

    H = nx.convert_node_labels_to_integers(G,first_label=0, ordering='default',label_attribute='old_label')

    options = {

    'node_color': 'lightblue',

    'node_size': 200,

    'width': 3,

    'arrowstyle': '-|>',

    'arrowsize': 10,

    }

    nx.draw_networkx(H,arrows=True,**options)

    return H

#Epigenomics Graph Generator

def Epigenomics(m):

    G = nx.DiGraph()

    num=1

    G.add_node(0);

    G.add_node(4*m+2-1)

    G.add_edge(4*m+2-1,4*m+3-1)

    G.add_edge(4*m+3-1,4*m+4-1)



    for k in range(1,m+1):

        prev = 0

        for i in range(1,5):

            G.add_edge(prev,num)

            prev = num

            num += m

        G.add_edge(prev,4*m+2)

        num = k+1

    return G

#Check category intersection
def CategoryConstraint(t_ij, a_ij):
    flag = False
    for i in C_a_i_j[a_ij-1]:
        for j in C_t_i_j[t_ij-1]:
            if i == j:
                flag = True
                #print(i,j,C_a_i_j[a_ij-1][i],C_t_i_j[t_ij-1][j])
                break
    return flag

def avail(a_i):
    #   let's say agent a_i is allocated to task t_k
    #   t_k = aloc[a_i]
    #   time taken for execution of t_k task to complete on agent a_i will give us the earliest time it is available for
    #   allocation and if the agent is not allocated then the array aloc[a_i] will have some appropriate value
    for x in aloc:
        if x == a_i:
            return (Time_ji[aloc[a_i]][a_i])
    return 0

def EST(t_j, a_i):
    if t_j == 0:
        return 0
    else :
        # req_Rank = Rank[t_j]-1
        # t_i=0
        # for i in range(1,numOfTasks+1):
        #    if Rank[i]==req_Rank:
        #     t_i = i
        #     break
        maxAFT = -np.inf
        predicessorsOfCurrentTask =  Graph.predecessors(t_j)
        for predicessor in predicessorsOfCurrentTask:
            if AFT(predicessor) > maxAFT:
                maxAFT = AFT(predicessor)
        return max(avail(a_i),maxAFT)

#def AST(t_j):
 #   return (EST(t_j, aloc[t_j]))

def AFT(t_j):
    return (EST(t_j,aloc[t_j]) + Time_ji[t_j][aloc[t_j]])

def Algo1(Graph):
    # Initializing a queue
    OpenN = Queue()
    # Adding elements to a queue
    t_exit = numOfTasks-1
    OpenN.put(t_exit)
    while not OpenN.empty() :
        t_j = OpenN.get()
        for a_i in agents:
            if t_j == t_exit :
                Zeta_ji[t_j][a_i] = Time_ji[t_j][a_i]
            else :
                Zeta_ji[t_j][a_i] = 0
                t_k = list(Graph.successors(t_j))
                for k in t_k :
                    X = np.inf
                    for r in agents :
                        X = min(X, Zeta_ji[k][r] + Time_ji[t_j][a_i])
                    Zeta_ji[t_j][a_i] = max(Zeta_ji[t_j][a_i], X)
        rank = 0
        for a_i in agents:
            rank += Zeta_ji[t_j][a_i]
        rank /= numOfAgents
        Rank[t_j] = rank

        maxSuccRank = 0
        t_k = list(Graph.successors(t_j))
        for k in t_k:
            maxSuccRank = max(maxSuccRank, Rank[k])

        if Rank[t_j] <= maxSuccRank :  
            for a_i in agents:
                Zeta_ji[t_j][a_i] = Zeta_ji[t_j][a_i] * (maxSuccRank + Delta)/Rank[t_j]
                Rank[t_j] = maxSuccRank + Delta

        predicessorsOfCurrentTask =  Graph.predecessors(t_j)
        for predicessor in predicessorsOfCurrentTask:
            OpenN.put(predicessor)                

def fill_mec(t_j, a_i, Graph):
    t_exit = numOfTasks-1
    if t_j == t_exit:
        MEC[t_j][a_i] = Time_ji[t_j][a_i]*cost[a_i]
    else :
        t_k = list(Graph.successors(t_j))
        Y = np.inf
        for k in t_k :
            for r in agents :
                Y = min(Y, MEC[k][r] + Time_ji[t_j][a_i] * cost[a_i])
        MEC[t_j][a_i] = Y


def schedule_task(Graph, numOfAgents, MEC, D):
    minlensched = True
    schedule = []
    # print("hsdjuksvl")
    while not t_s.empty():
        (a,t_j) = t_s.get()
        print(t_j)
        As =[]
        tempAloc = 0
        for a_i in agents:
            ESL[t_j][a_i]=D-(EST(t_j,a_i)+Zeta_ji[t_j][a_i])*1/tf
            if(ESL[t_j][a_i] > 0):
              As.append(a_i)
            # if not As.isEmpty():
        if len(As)!=0:
            ans = np.inf
            for a_i in As:
                if(ans > MEC[t_j][a_i]):
                    ans = MEC[t_j][a_i]
                    tempAloc = a_i
            aloc[t_j] = tempAloc
            minlensched = False
        else:
            ans = -np.inf
            for a_i in agents:
                if(ans < ESL[t_j][a_i]):
                    ans = ESL[t_j][a_i]
                    tempAloc = a_i
            aloc[t_j] = tempAloc
        schedule.append((t_j,aloc[t_j],EST(t_j,aloc[t_j])))
                       
    limit = AFT(numOfTasks-1)
    if limit <= D:
        totcost = 0
        for t_j in tasks: # doubt b/w agent n task
            totcost= totcost+Time_ji[t_j][aloc[t_j]]*cost[aloc[t_j]]
        return True, minlensched, totcost,schedule
                                                         
    return False, minlensched, int(0),schedule

def Algo2(Graph):
    # [Zeta,Rnk]=Algo1(Graph)
    for t_j in tasks[::-1]:
        for a_i in agents[::-1]:
            #print(str(t_j) +' '+ str(a_i))
            if MEC[t_j][a_i] == 0:
                fill_mec(t_j,a_i ,Graph)
   
    for i in tasks:
        t_s.put((-Rank[i],i))
        # print(t_s)
    valid_shed = False
    min_len_shed = False
    tf=1
   
    while (not valid_shed) and (not min_len_shed):
        [valid_shed, min_len_shed, total_cost,schedule] = schedule_task(Graph, numOfAgents, MEC, D)
        makespan = 0
        for j in agents:  #DOUBT B/W AGENT N TASK  
            makespan= makespan+Time_ji[j][aloc[j]]

        D_o = makespan-D

        if (D_o /D)*100 ==1:
            tf = tf - 0.01
        else:
            c=np.log2((D_o/D)*100)
            tf = tf - 0.01*c
    if makespan > D :
        for a in schedule:
            print(a)
           
        print('the given schedule exceeds the deadline with a deadline overshoot of '+str(D_o - D))
        print('total cost = ' + str(total_cost))
    else:
        for a in schedule:
            print(a)
        print('the given schedule does not exceed the deadline')
        print("Total Cost = " + str(total_cost))
       
#global variables
m = 0
numOfAgents = 3
numOfTasks = 0
#numOfTasks = 5
agents = [i for i in range(numOfAgents)]
tasks = []
Delta = 0.01
D = 150 # Deadline
tf = 1
t_s=PriorityQueue()
time = 0
C_a_i_j = []
for i in range(numOfAgents):
        num= np.random.randint(4,10)
        #col = np.random.randint(1,10,size=num+1)
        col = random.sample(range(1, 10), num)
        C_a_i_j.append(col)
C_t_i_j = []
Rank = []
MEC = []
aloc = []
ESL = [[]]
       

for i in range(1, 5):
    m = i*4
    numOfTasks = (m**2 + m -2)//2
    # numOfTasks =4*m+4
    print("num of tasks -- " + str(numOfTasks))
    tasks = []
    tasks = [i for i in range(numOfTasks)]
    Delta = 0.01
    D = 150 # Deadline
    tf = 1
    t_s=PriorityQueue()
    time = 0
    #Categories of tasks
    C_t_i_j.clear()
    for i in range(numOfTasks):
        num= np.random.randint(4,10)
        #col = np.random.randint(1,10,size=num+1)
        col = random.sample(range(1, 10), num)
        C_t_i_j.append(col)

    print("Categories of agents:")
    print(C_a_i_j)
    print("Categories of tasks:")
    print(C_t_i_j)


    # Time Matrix
    Time_ji = np.random.randint(1,10,size=(numOfTasks,numOfAgents))
    print("time:")
    print(Time_ji)

    # Zeta matrix
    rows, cols = (numOfTasks, numOfAgents)
    Zeta_ji = [[0 for i in range(cols)] for j in range(rows)]
    print("ZETA:")
    print(Zeta_ji)
   
    ESL = [[]]

    #Rank array
    Rank.clear()
    aloc.clear()
    #Rank = [0 for i in range(numOfTasks)]
    for i in range(numOfTasks):
        Rank.append(0)
        aloc.append(0)

    # Aloc
    #aloc = [0 for i in range(numOfTasks)]

    # ESL Matrix
    ESL = [[0 for i in range(cols)] for j in range(rows)]

    #cost array
    cost= np.random.randint(80,100, size=(numOfAgents))
    print("COST: ")
    print(cost)

    #MEC matrix
    rows, cols = (numOfTasks, numOfAgents)
    MEC.clear()
    MEC = [[0 for i in range(cols)] for j in range(rows)]
    print("initial mec:")
    print(MEC)
    Graph = Epigenomics(m)

    print("The successors of node 1 are:")

    print(list(Graph.successors(1)))
    print("Category Constraint : ")
    print(CategoryConstraint(1,1))
    Algo1(Graph)
    print("zeta ji : ")
    print(Zeta_ji)
    print("rank:")
    print(Rank)
    Algo2(Graph)
    print("MEC : ")
    print(MEC)
    print(" ")
    print(" ")
    print(" ")
    print(" ")
    print(" ")