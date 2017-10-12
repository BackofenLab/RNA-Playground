(function(root){
  // Creates namespaced functions within a closure.
  //
  // http://github.com/jweir/namespace
  // Copyright (c) 2011 John Weir john@famedriver.com
  //
  // License: Public Domain
  //
  // Usage (within a closure):
  //
  //   namepace(moduleName, function1, function2, ...);
  //
  // Example:
  //
  //   (function(){
  //     namespace("space", publicFn, anotherFn);
  //
  //     function privateFn(n) { return [n, "private"].join(" and "); }
  //     function publicFn(n)  { return privateFn(n); }
  //     function anotherFn(n) { return "another result"; }
  //   }());
  //
  //   space.publicFn("good"); /* => "good and private" */
  //
  // Accepts nested namespaces, even if the parent does not yet exist.
  //
  //   namespaces("foo.bar", func);
  //
  // If no string is given as the first argument, the functions will be scoped to the root.

  function namespace(moduleNamespace, functions) {
    var space, i;

    if (typeof arguments[0] != "string"){
      space = root;
      i = 0;
    } else {
      space = namespaceFor(arguments[0].split("."));
      i = 1;
    }

    for (var l = arguments.length; i < l; i++) {
      space[functionName(arguments[i])] = arguments[i];
    }
  }

  function functionName(fn) {
    return fn.name || fn.toString().match(/^\s*function\s+([^\s\(]+)/)[1];
  }

  function namespaceFor(parts) {
    var space = root;

    for (var i=0; i < parts.length; i++) {
      space = space[parts[i]] = space[parts[i]] || {};
    }
    return space;
  }

  namespace(namespace);
}(this)); // <- change this if you want to change the scope of namespace
