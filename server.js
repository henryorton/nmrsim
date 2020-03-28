const http = require('http');
const fs = require('fs');

function returnHomepage (response) {
  fs.readFile('dist/index.html', 
    function (error, content) {
      response.writeHead(200, {'Content-Type': 'text/html'});
      response.write(content);
      response.end();
    }
  );
}

function returnScripts (response) {
  fs.readFile('dist/main.js', 
    function (error, content) {
      response.writeHead(200, {'Content-Type': 'text/javascript'});
      response.write(content);
      response.end();
    }
  );
}

// Make server
http.createServer(
  function (request, response) {
    if (request.url == '/') {
      returnHomepage(response);
    }
    else if (request.url == '/main.js') {
      returnScripts(response);
    }
  }
).listen(8080);