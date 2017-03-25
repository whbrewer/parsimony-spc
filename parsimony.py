import ConfigParser
import math
import time

# A very simple representation for Nodes.
# Leaves are anything which is not a Node.
class Node(object):
  def __init__(self, left, right, data=None):
    self.left = left
    self.right = right

  def __repr__(self):
    return '(%s,%s)' % (self.left, self.right)

class Parsimony:
  mx = 0       # numerical label for indicating the current score node
  f = []       # set of possible base assignments
  l = []       # set of the number of changes
  topL = 0     # the left child node of the root
  topR = 0     # the right child node of the root
  flag = False # indicating the right part of the whole tree
  root = True  # indicate that the current executing point is just down of the root
  slen = 0
  species = []
  n = 0
  best = {'score':0,'id':0,'tree':(),'f_root':'','l_root':[]}
  outp = ''
  space = '&nbsp;'
  separator = '-'*80

  def __init__(self,species):
    self.species = species.split(',')
    self.slen = len(self.species[0])
    self.n = len(self.species)
    self.outp += 'PARSIMONY METHOD<br>'
    self.outp += self.separator + '<br>'
    self.f = [[ '' for i in range(self.slen) ] for i in range(2*self.n-1)]
    for i in range(self.n):
        self.outp += str(i) + ':' + self.space + str(self.species[i]) + self.space
    self.outp +='<br>' + self.separator + '<br>'
    for i in range(self.n):
        self.f[i] = list(self.species[i])
    self.l = [[ 0  for i in range(self.slen) ] for i in range(2*self.n-1)]

  def add_leaf(self,tree, label):
    '''Given a tree and a label, yields every possible augmentation of the
    tree by adding a new node with the label as a child "above" some
    existing Node or Leaf.'''
    yield Node(label, tree)
    if isinstance(tree, Node):
      for left in self.add_leaf(tree.left, label):
        yield Node(left, tree.right)
      for right in self.add_leaf(tree.right, label):
        yield Node(tree.left, right)

  def enum_unordered(self,labels):
     '''Given a list of labels, yield each rooted, unordered full binary tree
     with the specified labels.'''
     if len(labels) == 1:
       yield labels[0]
     else:
       for tree in self.enum_unordered(labels[1:]):
         for new_tree in self.add_leaf(tree, labels[0]):
           yield new_tree

  def computefl(self,f1,f2,l1,l2):
    '''compute F and L'''
    f3 = [ '' for i in range(self.slen) ]
    l3 = [ 0  for i in range(self.slen) ]
    # step through each item/character in list
    for i in range(len(f1)):
      intersect = self.intersection(f1[i],f2[i])
      if intersect:
        # all characters that get intersected get returned as
        # list of characters, so must join them together
        f3[i] = ''.join(intersect)
        l3[i] = l1[i] + l2[i]
      else:
        f3[i] = f1[i] + f2[i] #union
        l3[i] += l1[i] + l2[i] + 1
    return [f3,l3]

  def intersection(self,s1, s2):
    '''compute the intersection of two strings'''
    s1_chars = set(s1)
    s2_chars = set(s2)
    result = sorted(s1_chars.intersection(s2_chars))
    return result

  def parse(self,t):
    '''parse Newick tree recursively'''
    # base case(leave,leave)
    if type(t[0]) is int and type(t[1]) is int:
      a = t[0]; b = t[1]
      self.mx += 1

      if self.flag == True:   # indicates that the starting point
                            # is in the left part of the whole tree
          self.topL = self.mx # updates the topL value as mx
      else:
          self.topR = self.mx # updates the topR indicating the right part
                            # of the whole tree as mx
      #self.outp+='F='+str(self.mx)+self.space*10+str(t)+'<br>' # output pattern
      # compute F and L of the mx
      [self.f[self.mx],self.l[self.mx]] = self.computefl(self.f[a],self.f[b],\
                                                         self.l[a],self.l[b])
    elif type(t[0]) is tuple and type(t[1]) is int: # left part is subtree,
                                                    # right part is leave
      if self.root == True: # indicates the operations in the root case
         self.topR = t[1]    # topR is set as right leave of analyzing tree (in root case)
         self.root = False   # we don't consider the root case anymore next time
         self.flag = True    # we will consider the left subtree of the whole tree next time
      self.parse(t[0])      # reconsider the left subtree next time
      self.parse((t[1],self.mx)) # computes the F,L of mx+1 with left leave
                                 # and right leave
    elif type(t[0]) is int and type(t[1]) is tuple: # left part is leave,
                                                    # right part is subtree
      if self.root == True: # indicates the operations in the root case
         self.topL = t[0]   # topL is set as left leave of analyzing tree (in root case)
         self.root = False  # don't consider the root case anymore next time
         self.flag = False  # consider the right subtree of whole tree next time
      self.parse(t[1])      # reconsider the right subtree next time
      self.parse((t[0],self.mx)) # computes the F,L of mx+1 with left leave and right leave
    else: # tuple tuple case(left part is subtree, right part is subtree)
      if self.root == True:  # indicates the operations in the root case
        self.flag = True   # consider left subtree of whole tree next time
        self.root = False  # don't consider the root case anymore next time
        self.parse(t[0])   # reconsider the left subtree next time
        self.flag = False  # consider right subtree of whole tree
                           # only after analyzing left subtree of whole tree
        self.parse(t[1])
      else: #indicates the operations in the subtree cases, not the root case
        self.parse(t[0])   # reconsider the left subtree
        self.parse(t[1])   # reconsider the right subtree

  def max_trees(self,n):
    return math.factorial(2*n-3)/math.factorial(n-2)/2**(n-2)

  def result(self):
    start = time.clock()
    x = tuple([ str(i) for i in range(self.n) ])
    ii = 0
    self.best = { 'score': 2*self.n*self.slen }
    for tree in self.enum_unordered(x):
      self.mx = self.n-1 # index for f
      self.flag = False #initial value for flag
      self.topL = 0
      self.topR = 0
      self.root = True #initial value for indicating the root case
      ii += 1
      # convert from tree data structure to string
      s = str(tree)
      # convert from string to tuple
      t = eval(s)
      # parse tuple
      self.parse(t)
      # one last call to computefl for the final two
      if type(t[0]) is tuple and type(t[1]) is tuple:
        self.mx += 1 # mx must be increased in only case that the whole tree
                   # has a left subtree and a right subtree
      # e.g. if topR=8, topR must be changed into 7 for computing the final score
      if self.topR == 2*self.n-2: self.topR -= 1
      # e.g., if topL=8, topR must be changed into 7 for computing the final score
      if self.topL == 2*self.n-2: self.topL -= 1
      [self.f[self.mx],self.l[self.mx]] = self.computefl(self.f[self.topL],\
                                                         self.f[self.topR],\
                                                         self.l[self.topL],\
                                                         self.l[self.topR])
      score = sum(self.l[2*self.n-2])
      # store info for the best tree
      if score < self.best['score']:
        self.best['id'] = ii
        self.best['score'] = score
        self.best['tree'] = t
        self.best['f_root'] = ''.join(self.f[2*self.n-2])
        self.best['l_root'] = self.l[2*self.n-2]
      self.outp+=str(ii)+':'+self.space*5+str(s)+self.space*5+'score:'+str(score)+'<br>'

    self.outp += '-'*80 + '<br>'
    self.outp += 'BEST SCORING TREE:' + '<br>'
    self.outp += 'Id:'+self.space*5+str(self.best['id'])+'<br>'
    self.outp += 'Score:'+self.space*5+str(self.best['score'])+'<br>'
    self.outp += 'Tree:'+self.space*5+str(self.best['tree'])+'<br>'
    self.outp += 'F_root:'+self.space*5+str(self.best['f_root'])[0:self.best['l_root'].index(2)]+'['+str(self.best['f_root'])[self.best['l_root'].index(2):self.best['l_root'].index(2)+2]+']'+str(self.best['f_root'])[self.best['l_root'].index(2)+2:]+'<br>'
    self.outp+='L_root:'+self.space*5+str(self.best['l_root'])+'<br>'
    elapsed = time.clock() - start
    self.outp+= '-'*80 + '<br>'
    self.outp += 'Time elapsed:'+self.space*4+str(elapsed)[0:6]+'sec'

if __name__ == "__main__":
    # read parameters from input file
    params = dict()
    config = ConfigParser.ConfigParser()
    config.read('parsimony.ini')
    for section in config.sections():
        for option in config.options(section):
            params[option] = config.get(section, option)

    # parse out parameters from params dictionary
    species = params['species']
    par = Parsimony(species)
    par.result()
    print par.outp
