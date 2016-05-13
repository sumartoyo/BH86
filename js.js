var NODE_INTERNAL = 'nin';
var NODE_EXTERNAL = 'nex';
var THETA = 0.5;
var G = 0.0001;
var EPS_2 = 0.001;
var ETA = 0.1;

var Body = function(x, y, mass) {

    var self = this;

    self.x = x;
    self.y = y;
    self.mass = mass;
    self.vx = 0;
    self.vy = 0;

    self.computeForce = function(node) {

        var force = {
            x: 0,
            y: 0,
        };

        var dx = node.xmass-self.x;
        var dy = node.ymass-self.y;
        if (dx == 0 && dy == 0) {
            return force;
        }

        var r = Math.sqrt(Math.pow(dx, 2) + Math.pow(dy, 2));
        var strength = G*self.mass*node.mass /
            Math.sqrt(Math.pow(Math.pow(r, 2)+EPS_2, 3));

        if (node.body !== null) {
            force.x = strength*dx;
            force.y = strength*dy;
        } else {
            if (node.width/r < THETA) {
                force.x = strength*dx;
                force.y = strength*dy;
            } else {
                ['nw', 'ne', 'sw', 'se'].forEach(d => {
                    var child = node[d];
                    if (child !== null) {
                        var fChild = self.computeForce(child);
                        force.x += fChild.x;
                        force.y += fChild.y;
                    }
                });
            }
        }

        return force;
    };

    self.updatePos = function(force) {
        var ax = force.x/self.mass,
            ay = force.y/self.mass,
            dvx = ax*ETA,
            dvy = ay*ETA;
        self.x += (self.vx*ETA)+(0.5*dvx*ETA);
        self.y += (self.vy*ETA)+(0.5*dvy*ETA);
        self.vx += dvx;
        self.vy += dvy;
    };
};

var Node = function(xmin, ymin, width) {

    var self = this;

    self.type = NODE_EXTERNAL;
    self.body = null;
    self.mass = 0;
    self.xmass = 0;
    self.ymass = 0;

    self.xmin = xmin;
    self.ymin = ymin;
    self.width = width;

    self.nw = null;
    self.ne = null;
    self.sw = null;
    self.se = null;

    self.insert = function(body) {
        self.updateMass(body);
        if (self.type == NODE_INTERNAL) {
            self.insertRecursive(body);
        } else if (self.body !== null) {
            self.type = NODE_INTERNAL;
            self.branch();
            self.insertRecursive(self.body);
            self.insertRecursive(body);
            self.body = null;
        } else {
            self.body = body;
        }
    };

    self.insertRecursive = function(body) {
        var widthHalf = self.width/2;
        var insertTo = body.y < self.ymin+widthHalf ? 's' : 'n';
        insertTo += body.x < self.xmin+widthHalf ? 'w' : 'e';
        self[insertTo].insert(body);
    };

    self.updateMass = function(body) {
        var weightX = self.xmass*self.mass;
        var weightY = self.ymass*self.mass;
        weightX += body.x*body.mass;
        weightY += body.y*body.mass;
        self.mass += body.mass;
        self.xmass = weightX/self.mass;
        self.ymass = weightY/self.mass;
    };

    self.branch = function() {
        var widthHalf = self.width/2;
        self.nw = new Node(self.xmin, self.ymin+widthHalf, widthHalf);
        self.ne = new Node(self.xmin+widthHalf, self.ymin+widthHalf, widthHalf);
        self.sw = new Node(self.xmin, self.ymin, widthHalf);
        self.se = new Node(self.xmin+widthHalf, self.ymin, widthHalf);
    };

    self.eliminateEmpty = function() {
        ['nw', 'ne', 'sw', 'se'].forEach(d => {
            var child = self[d];
            if (child !== null) {
                if (child.type == NODE_EXTERNAL) {
                    if (child.body === null) {
                        self[d] = null;
                    }
                } else {
                    child.eliminateEmpty();
                }
            }
        });
    };

    self.print = function(tab, name) {
        console.log(tab+' '+name+' '+self.type+' '+(
            self.body !== null ?
                JSON.stringify(self.body) :
                'x='+self.xmass+' y='+self.ymass+' m='+self.mass
        ));
        ['nw', 'ne', 'sw', 'se'].forEach(d => {
            var child = self[d];
            if (child !== null) {
                child.print(tab+'\t', d);
            }
        });
    };
};

var Tree = function(bodies) {

    var xmin = null, xmax = null, ymin = null, ymax = null;
    bodies.forEach(b => {
        if (xmin === null || b.x < xmin) xmin = b.x;
        if (xmax === null || b.x > xmax) xmax = b.x;
        if (ymin === null || b.y < ymin) ymin = b.y;
        if (ymax === null || b.y > ymax) ymax = b.y;
    });

    var root = new Node(xmin, ymin, Math.max(xmax-xmin, ymax-ymin)+1);
    bodies.forEach(b => {
        root.insert(b);
    });
    root.eliminateEmpty();

    return root;
};

var Step = function(bodies) {

    var self = this;
    var isRun = true;

    self.bodies = bodies;

    self.run = function(epoch) {
        console.log('epoch', epoch);
        var tree = Tree(self.bodies);
        self.bodies.forEach(b => {
            var force = b.computeForce(tree);
            b.updatePos(force);
            console.log(b.mass+' '+b.x+' '+b.vx);
        });
        if (isRun) {
            setTimeout(function() {
                self.run(epoch+1);
            }, 100);
        }
    };

    self.stop = function() {
        isRun = false;
    };
};

var bodies = [
    // new Body(-3, -3, 3),
    new Body(-0.2, 0, 4),
    new Body(0.2, 0, 4),
];
// var root = Tree(bodies);
// root.print('', 'root');
var step = new Step(bodies);

function r() {
    window.location.reload();
}

function run() {
    step.run(1);
}
