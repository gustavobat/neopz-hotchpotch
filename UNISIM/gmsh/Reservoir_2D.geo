cl__1 = 1e+22;
Point(1) = {-379.1033, -125.26547, 0, cl__1};
Point(2) = {769.9057, -221.73367, 0, cl__1};
Point(3) = {-379.1033, -225.26547, 0, cl__1};
Point(4) = {769.9057, -321.73367, 0, cl__1};
Point(5) = {-379.1033, -125.26547, 0, cl__1};
Point(6) = {-379.1033, -225.26547, 0, cl__1};
Point(7) = {769.9057, -221.73367, 0, cl__1};
Point(8) = {769.9057, -321.73367, 0, cl__1};
p1 = newp;
Point(p1 + 1) = {-315.7472206646397, -99.45325798383324, 0};
Point(p1 + 2) = {-268.7720603736423, -50.46644353771339, 0};
Point(p1 + 3) = {-237.6034270446602, 9.489625031556697, 0};
Point(p1 + 4) = {-190.2198094696481, 40.80690719931638, 0};
Point(p1 + 5) = {-120.9858283296585, 39.86528192416553, 0};
Point(p1 + 6) = {-76.26573323287643, 84.04718188882102, 0};
Point(p1 + 7) = {-23.2616535304749, 121.9087676382769, 0};
Point(p1 + 8) = {38.56305545179947, 123.8801210713907, 0};
Point(p1 + 9) = {99.68311319254995, 93.21466364041368, 0};
Point(p1 + 10) = {155.7199063135849, 55.97507191806057, 0};
Point(p1 + 11) = {220.7635187902494, 42.96567763464536, 0};
Point(p1 + 12) = {270.785932489392, -6.768071589551177, 0};
Point(p1 + 13) = {329.4806621823897, -41.41491596324465, 0};
Point(p1 + 14) = {391.7374492423589, -69.80768884257547, 0};
Point(p1 + 15) = {450.8936931968585, -103.752069775442, 0};
Point(p1 + 16) = {510.7002203654192, -136.3712266373657, 0};
Point(p1 + 17) = {576.6039340270219, -155.1230600147117, 0};
Point(p1 + 18) = {640.4222305054384, -180.2144755433361, 0};
Point(p1 + 19) = {707.1847965251944, -195.0635847322544, 0};

p2 = newp;
Point(p2 + 1) = {-315.7472206646397, -199.4532579459426, 0};
Point(p2 + 2) = {-268.7720603736423, -150.4664433665869, 0};
Point(p2 + 3) = {-237.6034270446602, -90.51037500218166, 0};
Point(p2 + 4) = {-190.2198094696481, -59.19309280068363, 0};
Point(p2 + 5) = {-120.9858283296585, -60.13471807583447, 0};
Point(p2 + 6) = {-76.26573323287643, -15.95281827271342, 0};
Point(p2 + 7) = {-23.2616535304749, 21.9087674791631, 0};
Point(p2 + 8) = {38.56305545179947, 23.88012070369869, 0};
Point(p2 + 9) = {99.68311319254995, -6.785336484251856, 0};
Point(p2 + 10) = {155.7199063135849, -44.02492808187778, 0};
Point(p2 + 11) = {220.7635187902494, -57.03432234700983, 0};
Point(p2 + 12) = {270.785932489392, -106.7680714818664, 0};
Point(p2 + 13) = {329.4806621823897, -141.4149158144327, 0};
Point(p2 + 14) = {391.7374492423589, -169.8076888038037, 0};
Point(p2 + 15) = {450.8936931968585, -203.7520697754323, 0};
Point(p2 + 16) = {510.7002203654192, -236.3712266373657, 0};
Point(p2 + 17) = {576.6039340270219, -255.1230600147117, 0};
Point(p2 + 18) = {640.4222305054384, -280.2144755433361, 0};
Point(p2 + 19) = {707.1847965251944, -295.0635847322544, 0};

If (xzQ == 1)

// Apply rotation
rotate_p[] = Rotate { { 1, 0,0}, {0, 0, 0}, -Pi/2.0 } {
Point {1, p1 + 1, p1 + 2, p1 + 3, p1 + 4, p1 + 5, p1 + 6, p1 + 7, p1 + 8, p1 + 9, p1 + 10, p1 + 11, p1 + 12, p1 + 13, p1 + 14, p1 + 15, p1 + 16, p1 + 17, p1 + 18, p1 + 19, 2, 3, p2 + 1, p2 + 2, p2 + 3, p2 + 4, p2 + 5, p2 + 6, p2 + 7, p2 + 8, p2 + 9, p2 + 10, p2 + 11, p2 + 12, p2 + 13, p2 + 14, p2 + 15, p2 + 16, p2 + 17, p2 + 18, p2 + 19, 4};
};

};

EndIf

Spline(1) = {1, p1 + 1, p1 + 2, p1 + 3, p1 + 4, p1 + 5, p1 + 6, p1 + 7, p1 + 8, p1 + 9, p1 + 10, p1 + 11, p1 + 12, p1 + 13, p1 + 14, p1 + 15, p1 + 16, p1 + 17, p1 + 18, p1 + 19, 2};

Spline(2) = {3, p2 + 1, p2 + 2, p2 + 3, p2 + 4, p2 + 5, p2 + 6, p2 + 7, p2 + 8, p2 + 9, p2 + 10, p2 + 11, p2 + 12, p2 + 13, p2 + 14, p2 + 15, p2 + 16, p2 + 17, p2 + 18, p2 + 19, 4};
Line(3) = {3, 1};
Line(4) = {4, 2};

