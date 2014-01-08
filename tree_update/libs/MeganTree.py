#!/usr/bin/python

try:
    #from python_modules.LCAComputation import *
    from LCAComputation import *
    from os import sys
except:
    print """ Could not load some user defined module functions"""
    sys.exit(3)


class MeganTree:
    begin_pattern = re.compile("#")

    child_to_parent = {}
    parent_to_child = {}
    id_to_name = {}
    taxid_to_ptaxid = {}
    children_at_leaf = {}
    lca = None
    output = None
    def __init__(self, lca):
       self.lca = lca
       # TODO: Bring LCA Init to this init
    
    ### Getter functions
    def getChildToParentMap(self):
       return self.child_to_parent

    def getParentToChildrenMap(self):
       if not self.parent_to_child:
         self.computeParentToChildrenMap()
       return self.parent_to_child
    
    # function to build tree for all parent and child pairs
    def build_tree(self):
        print "In build_tree()"
        # for each child_id in taxid_to_ptaxid
        for t_id in self.lca.taxid_to_ptaxid:
             self.insertTaxon(t_id)
        print "Finished build_tree()"

    # function to insert taxon to the top-down (partent-to-child) mapping
    def insertTaxon(self, oid):
       id = oid 
       # climb the chain of parents until you are at the root node (id=1)
       while id != '1':
          pid = self.lca.getParentTaxId(id) # get the parent
          if pid is None:
             return
          if not id in self.child_to_parent:
             self.child_to_parent[id] = [pid, 0]
          self.child_to_parent[id][1] += 1
          id = pid

       if not oid in self.children_at_leaf:
         self.children_at_leaf[oid] = 0
       self.children_at_leaf[oid] +=1

    def computeParentToChildrenMap(self):
        print "In computeParentToChildrenMap:"
        for id, pid in self.child_to_parent.iteritems(): 
            if not pid[0] in self.parent_to_child:
                 self.parent_to_child[pid[0]] = [ [],  pid[1] ]
            self.parent_to_child[pid[0]][0].append(id)
        print "finished"

    def printTree(self, id):
        if not self.parent_to_child:
          self.computeParentToChildrenMap()

        self.output = ""
        self.createTree(id)
        return self.output
        
    def createTree(self, id):
        if not id in self.parent_to_child: 
           # leaf node
           if id=='1':
              print self.parent_to_child
              sys.exit(0)
           self.output +=  str(self.lca.translateIdToName(id))
           # self.output +=  id + ":" +  str(0)
           return 
        
        children =  self.parent_to_child[id][0] 
        count =  self.parent_to_child[id][1] 
        
        # self.output += "("
        i = 0
        self.output += "("
        for child in children:
            if i > 0 :
               self.output += ","
            self.createTree(child)
            i+=1
            if (i % 10000) == 0:
                print "i = " + str(i)
        self.output += ")"
        # self.output += ")" +  str(self.lca.translateIdToName(id)) + "}" +  str(count)
        self.output += str(self.lca.translateIdToName(id))
        if id == "1":
            self.output += ";"
        ##self.output += ")" +  id + ":" +  str(count)
    
    # driver function for recurisive drawNewickTree
    def printNewickTree(self, id):
        # check to ensure that parent to child mapping present
        if not self.parent_to_child:
            self.computeParentToChildrenMap()
        self.output = "" # reset output variable
        self.createNewickTree(id)
        return self.output # return tree string
    
    # prints a root node (id=1) prints a tree in the Newick format
    # TODO: could be generalized to print subtree from any starting node
    # TODO: this function is very slow, should be some way to speed it up
    def createNewickTree(self, id):
        # determine if at leaf
        if not id in self.parent_to_child:
            # at leaf print id
            self.output += str(id)
            return
        # get children of current node
        children = self.parent_to_child[id][0]
        
        # print children recursively
        i = 0 # children counter
        self.output += "("
        for child in children:
            # after first child print commas
            if i > 0 :
               self.output += ","
            self.createNewickTree(child) # recursive call
            i += 1
        self.output += ")" # close subtree
        self.output += str(id) # print internal nodes
        # if you are the root node finish with semicolon
        if id == "1":
            self.output += ";"
        return
