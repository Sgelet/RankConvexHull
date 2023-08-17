//
// Created by etoga on 1/16/23.
//

#ifndef CONVEX_HULL_H
#define CONVEX_HULL_H

#include <cstddef>
#include <algorithm>
#include <iostream>
#include <tuple>
#include <vector>
#include <deque>

template<class T>
class CHTree {
    struct Node;

    struct Point{
        double x,y;
    };

    struct Line{
        double slope,offset;

        double Eval(const double& x){ return slope*x+offset; }
    };

    struct Leaf {
        Leaf *pre = nullptr, *nxt = nullptr;
        T val;

        explicit Leaf(T v) : val(v) {};

        Leaf(T v, Leaf *l, const bool insertPre) : val(v) {
            if (!l) return;

            if (insertPre) {
                auto t = l->pre;
                l->pre = this;
                nxt = l;
                pre = t;
                if (t) t->nxt = this;
            } else {
                auto t = l->nxt;
                l->nxt = this;
                nxt = t;
                pre = l;
                if (t) t->pre = this;
            }
        }

        ~Leaf(){
            if(pre) pre->nxt = nxt;
            if(nxt) nxt->pre = pre;
        }
    };


    struct Bridge{
        T l,r;
        size_t ldiff,rdiff,tdiff;
    };

    struct Node {
        Node *left = nullptr, *right = nullptr, *par = nullptr;
        Leaf *min,*max;
        Bridge uBridge;
        Bridge lBridge;
        size_t size = 1, rank = 1;
        int balance = 0;

        Node() = default;

        explicit Node(Leaf* v) : min(v),max(v) {};

        ~Node(){
            if(!left && !right) delete min;
            if(left) delete left;
            if(right) delete right;
        }
    };


public:
    void Insert(T key){
        if(!root){
            Leaf* leaf = new Leaf(key);
            Node* node = new Node(leaf);
            InitBridge(node);
            root = node;
            return;
        }

        auto leaf = Find(key);

        if(leaf->min->val == key) return;

        Leaf* old_leaf = leaf->min;

        Node* left;
        Node* right;

        if(key<old_leaf->val){
            left = new Node(new Leaf(key, old_leaf, true));
            right = new Node(old_leaf);
        } else {
            left = new Node(old_leaf);
            right = new Node(new Leaf(key, old_leaf, false));
        }


        MakeLeftChild(left,leaf);
        MakeRightChild(right,leaf);

        InitBridge(left);
        InitBridge(right);

        left->rank = leaf->rank;
        leaf->rank++;
        right->rank = leaf->rank;

        Update(leaf,true, key);
    }

    void Remove(T key){
        if(!root) return;

        auto leaf = Find(key);

        if(leaf->min->val != key) return;

        // Remove current and link
        if(!leaf->par){
            delete leaf;
            root = nullptr;
            return;
        }

        // Replace parent with sibling
        Node *sibling;
        if(IsLeftChild(leaf)){
            sibling = leaf->par->right;
            leaf->par->right = nullptr;
        }
        else {
            sibling = leaf->par->left;
            leaf->par->left = nullptr;
        }

        Node *grandparent = leaf->par->par;
        sibling->par = grandparent;
        if(!grandparent){
            delete leaf->par;
            root = sibling;
            return;
        }
        if(IsLeftChild(leaf->par)) MakeLeftChild(sibling,grandparent);
        else MakeRightChild(sibling,grandparent);
        delete leaf->par;

        // Rebalance
        Update(sibling,false, key);
    }

    bool Covers(double x, double y){
        // Find bridge that covers upper and check if we lie under
        auto bridge_node = FindUCovering(x);
        if(!bridge_node) return false;
        auto[eu,u] = ExtractUBridge(bridge_node);
        if(eu.Eval(x) < y) return false;

        // Find bridge that covers lower and check if we lie above
        bridge_node = FindLCovering(x);
        if(!bridge_node) return false;
        auto[el, l] = ExtractLBridge(bridge_node);
        if(el.Eval(x) > y) return false;

        return true;
    }

    void PrettyPrinter(){
        PrintHelper(root,"",true);
    }

    void Verify(){
        Checker(root);
    }

    std::vector<std::pair<int,int>> UpperPoints(){
        if(!root) return {};
        std::deque<std::pair<int,int>> dq;
        root->rank = size(root->left) + 1;
        dq.push_front({root->rank-root->uBridge.ldiff, root->uBridge.l});
        dq.push_back({root->rank+root->uBridge.rdiff, root->uBridge.r});
        while(dq.front().first!=1){
            auto pre = FindFirstU(dq.front().second,true);
            dq.push_front({pre->rank-pre->uBridge.ldiff, pre->uBridge.l});
        }
        while(dq.back().first!=root->size){
            auto succ = FindFirstU(dq.back().second,false);
            dq.push_back({succ->rank+succ->uBridge.rdiff, succ->uBridge.r});
        }
        std::vector<std::pair<int,int>> res;
        for(auto e: dq){
            res.push_back(e);
        }
        return res;
    }
    std::vector<std::pair<int,int>> LowerPoints(){
        if(!root) return {};
        std::deque<std::pair<int,int>> dq;
        root->rank = size(root->left) + 1;
        dq.push_front({root->rank-root->lBridge.ldiff, root->lBridge.l});
        dq.push_back({root->rank+root->lBridge.rdiff, root->lBridge.r});
        while(dq.front().first!=1){
            auto pre = FindFirstL(dq.front().second,true);
            dq.push_front({pre->rank-pre->lBridge.ldiff, pre->lBridge.l});
        }
        while(dq.back().first!=root->size){
            auto succ = FindFirstL(dq.back().second,false);
            dq.push_back({succ->rank+succ->lBridge.rdiff, succ->lBridge.r});
        }
        std::vector<std::pair<int,int>> res;
        for(auto e: dq){
            res.push_back(e);
        }
        return res;
    }

    ~CHTree(){
       if(root) delete root;
    }

protected:
    Node *root = nullptr;

    int Checker(Node *x){
        if(!x) return 0;

        auto checkL = Checker(x->left);
        auto checkR = Checker(x->right);

        // Checker(right) - Checker(left) == balance factor
        if(x->balance != checkR - checkL){
            std::cout << "INCORRECT BALANCE" << std::endl;
        }

        // Balance factor -1,0,1
        if(x->balance < -1 || x->balance > 1){
            std::cout << "TREE NOT BALANCED" << std::endl;
        }
        // Check sizes

        if(x->size != size(x->left)+ size(x->right) + IsLeaf(x)){
            std::cout << "INCORRECT SIZE" << std::endl;
        }
        // Check min,max

        return std::max(checkL,checkR) + 1;
    }

    void PrintHelper(Node* t, std::string indent, bool last){
        if(!t) return;
        std::cout << indent;
        if(last){
            std::cout << "R----";
            indent += "     ";
        } else {
            std::cout << "L----";
            indent += "|    ";
        }

        std::cout << "[" << t->uBridge.l << "," << t->uBridge.r << "][" << t->uBridge.ldiff << "," << t->uBridge.rdiff << "]( BF = " << t->balance << ", size = " << t->size << " )" << std::endl;

        PrintHelper(t->left, indent, false);
        PrintHelper(t->right, indent, true);
    }

    Node* Find(T key){
        auto current = root;
        current->rank = size(current->left) + 1;
        while(!IsLeaf(current)){
            if(current->right->min->val <= key) current = StepRight(current);
            else current = StepLeft(current);
        }
        return current;
    }

    Node* FindUCovering(double x){
        auto current = root;
        if(!root) return nullptr;
        root->rank = size(root->left) + 1;
        while(current && !(current->rank-current->uBridge.ldiff <= x && x <= current->rank+current->uBridge.rdiff)){
            if(x < current->rank - current->uBridge.ldiff ) current = StepLeft(current);
            else if (x > current->rank + current->uBridge.rdiff) current = StepRight(current);
            else return nullptr;
        }
        return current;
    }
    Node* FindLCovering(double x){
        auto current = root;
        if(!root) return nullptr;
        root->rank = size(root->left) + 1;
        while(current && !(current->rank-current->lBridge.ldiff <= x && x <= current->rank+current->lBridge.rdiff)){
            if(x < current->rank - current->lBridge.ldiff) current = StepLeft(current);
            else if (x > current->rank + current->lBridge.rdiff) current = StepRight(current);
            else return nullptr;
        }
        return current;
    }

    Node* FindFirstU(T key, bool left){
        auto current = root;
        current->rank = size(current->left) + 1;
        while((!left && current->uBridge.l != key) || (left && current->uBridge.r != key)){
            if(current->right->min->val <= key) current = StepRight(current);
            else current = StepLeft(current);
        }
        return current;
    }
    Node* FindFirstL(T key, bool left){
        auto current = root;
        size_t rank = size(current->left) + 1;
        while((!left && current->lBridge.l != key) || (left && current->lBridge.r != key)){
            if(current->right->min->val <= key) current = StepRight(current);
            else current = StepLeft(current);
        }
        return current;
    }
    inline
    size_t size(Node *x){
        return x ? x->size : 0;
    }

    inline
    bool IsLeftChild(Node *x){
        return x == x->par->left;
    }

    inline
    bool IsLeaf(Node *x){
        return x && !x->left && !x->right;
    }

    inline
    void MakeLeftChild(Node *x, Node *y){
        x->par = y;
        if(y) y->left = x;
    }

    inline
    void MakeRightChild(Node* x, Node *y){
        x->par = y;
        if(y) y->right = x;
    }

    inline
    void UpdateMinMax(Node *x){
        x->min = x->left->min;
        x->max = x->right->max;
    }

    void RotateRight(Node* x){
        Node *y = StepLeft(x);

        MakeLeftChild(y->right,x);

        if(!x->par){
            root = y;
            y->par = nullptr;
        }
        else if(IsLeftChild(x)) MakeLeftChild(y,x->par);
        else MakeRightChild(y,x->par);

        MakeRightChild(x,y);

        y->size += size(x->right);
        x->size -= size(y->left);

        UpdateMinMax(x);
        UpdateMinMax(y);
        ComputeBridge(x);
        ComputeBridge(y);

        x->balance += 1 - std::min(0, y->balance);
        y->balance += 1 + std::max(0, x->balance);
    }

    void RotateLeft(Node* x){
        Node *y = StepRight(x);

        MakeRightChild(y->left,x);

        if(!x->par){
            root = y;
            y->par = nullptr;
        }
        else if(IsLeftChild(x)) MakeLeftChild(y,x->par);
        else MakeRightChild(y,x->par);

        MakeLeftChild(x,y);

        y->size += size(x->left);
        x->size -= size(y->right);

        UpdateMinMax(x);
        UpdateMinMax(y);
        ComputeBridge(x);
        ComputeBridge(y);

        x->balance += - 1 - std::max(0, y->balance);
        y->balance += - 1 + std::min(0, x->balance);
    }

    void Rebalance(Node *x){
        if(x->balance > 0){
            if(x->right->balance < 0) RotateRight(x->right);
            RotateLeft(x);
        } else if (x->balance < 0){
            if(x->left->balance > 0) RotateLeft(x->left);
            RotateRight(x);
        }
    }

    void Update(Node *x, const bool insertion, const T val){
        bool balance = true;

        // Traverse to root
        for(auto current = x; current; current = current->par){
            if(!IsLeaf(current)) {
                current->size = size(current->left) + size(current->right);
                UpdateMinMax(current);
                ComputeBridge(current);
            }
            if(!balance){
                continue;
            }
            if(current->balance < -1 || current->balance > 1){
                Rebalance(current);
                current = current->par;
                if(insertion || (current->balance == -1 || current->balance == 1)){
                    balance = false;
                    continue;
                }
            }
            if(current->par){
                current->par->balance += 1 - 2*(insertion == IsLeftChild(current));
                if ((insertion && current->par->balance == 0) || (!insertion && (current->par->balance == -1 || current->par->balance == 1))) balance = false;
            }
        }
    }



    bool ContainTest(Bridge* b, T val, size_t relative_rank){
        double slope = (b->r - b->l)/b->totaldiff;
        double offset = b->l - slope;
        double fl = slope*relative_rank + offset;
        return val <= fl;
    }

    Bridge FindUpperBridge(Node* v){
        Node* x = StepLeft(v);
        Node* y = StepRight(v);
        Line e_l, e_r, lr;
        Point l,r;
        bool undecided;
        while (!(IsLeaf(x) && IsLeaf(y))) {
            undecided = true;
            std::tie(e_l, l) = ExtractUBridge(x);
            std::tie(e_r, r) = ExtractUBridge(y);
            lr = LineBetweenPoints(l,r);
            if (e_l.slope <= lr.slope && !IsLeaf(x)){
                x = StepLeft(x); undecided = false;
            }
            //else if (IsLeaf(y)){
            //    x  = StepRight(x);
            //}
            if (lr.slope <= e_r.slope && !IsLeaf(y)) {
                y = StepRight(y); undecided = false;
            }
            //else if (IsLeaf(x) && undecided){
            //    y = StepLeft(y);
            //}
            if (undecided) { // Case i
                double m = x->rank + size(x->right) - 0.5;
                auto a = e_l.Eval(m);
                auto b = e_r.Eval(m);
                if (!IsLeaf(x) && a >= b || IsLeaf(y)) {
                    x = StepRight(x);
                } else {
                    y = StepLeft(y);
                }
            }
        }
        return {x->min->val,y->max->val,v->rank - x->rank, y->rank - v->rank, y->rank - x->rank};
    }
    // TODO: Merge this into one
    Bridge FindLowerBridge(Node* v){
        Node* x = StepLeft(v);
        Node* y = StepRight(v);
        Line e_l, e_r, lr;
        Point l,r;
        bool undecided;
        while (!(IsLeaf(x) && IsLeaf(y))) {
            undecided = true;
            std::tie(e_l, l) = ExtractLBridge(x);
            std::tie(e_r, r) = ExtractLBridge(y);
            lr = LineBetweenPoints(l,r);
            if (e_l.slope >= lr.slope && !IsLeaf(x)){
                x = StepLeft(x); undecided = false;
            }
            //else if (IsLeaf(y)){
            //    x  = StepRight(x);
            //}
            if (lr.slope >= e_r.slope && !IsLeaf(y)) {
                y = StepRight(y); undecided = false;
            }
            //else if (IsLeaf(x) && undecided){
            //    y = StepLeft(y);
            //}
            if (undecided) { // Case i
                double m = x->rank + size(x->right) - 0.5;
                auto a = e_l.Eval(m);
                auto b = e_r.Eval(m);
                if (!IsLeaf(x) && a <= b || IsLeaf(y)) {
                    x = StepRight(x);
                } else {
                    y = StepLeft(y);
                }
            }
        }
        return {x->min->val,y->max->val,v->rank - x->rank, y->rank - v->rank, y->rank - x->rank};
    }

    Node* StepLeft(Node* v){
        if(!v->left) return nullptr;
        auto x = v->left;
        x->rank = (v->rank - size(x->right) - IsLeaf(x));
        return x;
    }

    Node* StepRight(Node* v){
        if(!v->right) return nullptr;
        auto x = v->right;
        x->rank = v->rank + size(x->left);
        return x;
    }

    Line LineBetweenPoints(const Point& v0, const Point& v1){
        Line l;
        l.slope = (v1.y - v0.y)/(v1.x - v0.x);
        l.offset = v0.y - l.slope * v0.x;
        return l;
    }


    std::pair<Line,Point> ExtractUBridge(Node* v){
        if(IsLeaf(v)) return {{0,0},{v->rank,v->min->val}};
        Line el = LineBetweenPoints({v->rank-v->uBridge.ldiff, v->uBridge.l}, {v->rank + v->uBridge.rdiff, v->uBridge.r});
        auto x = v->rank - v->uBridge.ldiff + (double)v->uBridge.tdiff / 2;
        return {el,{x,el.Eval(x)}};
    }

    //TODO: Merge these?
    std::pair<Line,Point> ExtractLBridge(Node* v){
        if(IsLeaf(v)) return {{0,0},{v->rank,v->min->val}};
        Line el = LineBetweenPoints({v->rank-v->lBridge.ldiff, v->lBridge.l}, {v->rank + v->lBridge.rdiff, v->lBridge.r});
        auto x = v->rank - v->lBridge.ldiff + (double)v->lBridge.tdiff / 2;
        return {el,{x,el.Eval(x)}};
    }
    void InitBridge(Node* v){
        v->uBridge = {v->min->val, v->max->val, 0, 0, 0};
        v->lBridge = {v->min->val, v->max->val, 0, 0, 0};
    }

    void ComputeBridge(Node* v){
        v->uBridge = FindUpperBridge(v);
        v->lBridge = FindLowerBridge(v);
    }
};
#endif //CONVEX_HULL_H
