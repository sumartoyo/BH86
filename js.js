var NIN = 99;
var NEX = -99;

var Body = function(x, y, mass) {

    var self = this;

    self.x = x;
    self.y = y;
    self.mass = mass;
};

var Node = function(startX, startY, width) {

    var self = this;

    self.type = NEX;
    self.bodies = [];
    self.mass = 0;
    self.centerX = -1;
    self.centerY = -1;

    self.startX = startX;
    self.startY = startY;
    self.width = width;

    self.nw = null;
    self.ne = null;
    self.sw = null;
    self.se = null;

    self.insert = function(body) {
        if (self.bodies.length == 0) {
            self.bodies.push(body);
            self.updateMass(body);
        } else if (self.type == NIN) {
            self.bodies.push(body);
            self.updateMass(body);
            self.insertRecursive(body);
        } else if (self.type == NEX) {
            self.bodies.push(body);
            self.updateMass(body);
            self.type = NIN;
            self.branch();
            self.bodies.forEach(b => {
                self.insertRecursive(b);
            });
        }
    };

    self.insertRecursive = function(body) {
        var widthHalf = self.width/2;
        var insertTo = self.startY <= body.y && body.y < self.startY+widthHalf ? 'n' : 's';
        insertTo += self.startX <= body.x && body.x < self.startX+widthHalf ? 'w' : 'e';
        self[insertTo].insert(body);
    };

    self.updateMass = function(body) {
        self.mass += body.mass;
        var centerX = 0;
        var centerY = 0;
        self.bodies.forEach(b => {
            centerX += b.x * b.mass;
            centerY += b.y * b.mass;
        });
        self.centerX = centerX/self.mass;
        self.centerY = centerY/self.mass;
    };

    self.branch = function() {
        var widthHalf = self.width/2;
        self.nw = new Node(self.startX, self.startY, widthHalf);
        self.ne = new Node(self.startX+widthHalf, self.startY, widthHalf);
        self.sw = new Node(self.startX, self.startY+widthHalf, widthHalf);
        self.se = new Node(self.startX+widthHalf, self.startY+widthHalf, widthHalf);
    };
};

var root = new Node(0, 0, 8);
