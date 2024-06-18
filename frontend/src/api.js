export var backend_url = "";
if(process.env.NODE_ENV === 'production') {
    var backend_url = process.env.REACT_APP_API_URL;
}


export function getHDNeighbors(selectedPoints, kNeighbors, metric) {
    return new Promise((resolve, reject) => {
        if (selectedPoints.length > 0) {
            fetch(`${backend_url}/backend/computeHDNeighbors`, {
                method: "POST",
                body: JSON.stringify({
                    k: kNeighbors,
                    points: selectedPoints,
                    hd_metric: metric,
                }),
                headers: {
                    "Content-type": "application/json; charset=UTF-8"
                }
            })
                .then(response => {
                    if (!response.ok) {
                        console.log(`/backend/computeHDNeighbors got HTTP error: Status ${response.status}, Text ${response.statusText}`);
                        resolve(response);
                    } else {
                        response.json().then(data => {
                            let newOpacities = data["binary"]
                            resolve(newOpacities);
                        })
                    }
                })
            }
        else {
            resolve([]);
        }
    })
}